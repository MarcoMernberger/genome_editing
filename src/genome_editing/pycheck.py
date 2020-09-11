#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pysam.py: Contains ...."""

from typing import List, Dict, Tuple
from mbf_align.lanes import Sample
from mbf_genomes import GenomeBase
from pysam import AlignedSegment
from pandas import DataFrame
from pypipegraph import Job
from pathlib import Path
import pandas as pd
import pypipegraph as ppg
import pysam
import numpy as np
import mbf_genomes
import mbf_genomics
import collections
import mbf_r
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as numpy2ri

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class PysamCheck:
    """The idea is to use pysam to check the cigar strings of alinged reads to
    extract the editing event"""

    def __init__(
        self,
        primers_or_range_by_amplicon_id: Dict[str, List[str]],
        event_window_by_amplicon: Dict[str, List[int]],
        genome: GenomeBase,
        paired: bool = True,
        mapping_qality_filter: int = 20,
        single_read_sufficient: bool = False,
        types_to_consider: List[str] = None,
    ):
        self.name = "PysamCheck"
        self.primers_or_range_by_amplicon = primers_or_range_by_amplicon_id  # events within primer ranges are not reported
        self.event_window_by_amplicon = (
            event_window_by_amplicon  # reads need to span the event window
        )
        self.genome = genome
        self.paired = paired
        self.discard_if_pairs_differ = True
        self.single_read_sufficient = single_read_sufficient
        if types_to_consider is None:
            self.types_to_consider = ["insertion", "deletion", "mismatch"]

    def get_dependencies(self) -> List[Job]:
        deps = [
            ppg.FunctionInvariant("map_pairwise_space_to_reference_space", self.map_pairwise_space_to_reference_space),
            ppg.FunctionInvariant("evaluate_cigar", self.evaluate_cigar),
            ppg.FunctionInvariant("tally_reads", self.tally_reads),
            ppg.FunctionInvariant("map_events_to_reference_space", self.map_events_to_reference_space),
            ppg.FunctionInvariant("filter_events", self.filter_events),
            ppg.FunctionInvariant("extract_events", self.extract_events),
            ppg.FunctionInvariant("__check_reads_spanning_windows", self.__check_reads_spanning_windows),
            ppg.FunctionInvariant("__set_primer_range", self.__set_primer_range),
            ppg.FunctionInvariant("_prepare_reference_lookup", self._prepare_reference_lookup),
            self.genome.download_genome()
        ]
        return deps

    def _prepare_reference_lookup(self):
        """put the reference sequences in a dict for quick lookup"""

        def __do_lookup():
            lookup: Dict[str, str] = {}
            for name, length in self.genome.get_chromosome_lengths().items():
                lookup[name] = self.genome.get_genome_sequence(name, 0, length)
            return lookup

        job = ppg.AttributeLoadingJob(f"{self.name}__lookup", self, "_ref_lookup", __do_lookup).depends_on(self.get_dependencies())
        return job

    def __set_primer_range(self):
        """set the range in which we report events"""

        def __set_range_from_primers():
            if isinstance(next(iter(self.primers_or_range_by_amplicon.values()))[0], str):
                ranges: Dict[str, Tuple(int)] = {}
                for amplicon in self.primers_or_range_by_amplicon:
                    primer_forward, primer_reverse = self.primers_or_range_by_amplicon[
                        amplicon
                    ]
                    seq = self._ref_lookup[amplicon]
                    i1 = seq.find(primer_forward) + len(primer_forward)
                    i2 = seq.find(primer_reverse) + len(primer_reverse)
                    ranges[amplicon] = (i1, i2)
                return ranges
            else:
                return self.primers_or_range_by_amplicon

        job = ppg.AttributeLoadingJob(f"{self.name}__ranges_by_amplicon", self, "ranges_by_amplicon", __set_range_from_primers)
        job.depends_on(self.get_dependencies())
        job.depends_on(self._prepare_reference_lookup())
        return job

    def __check_reads_spanning_windows(self, reads: List[AlignedSegment]):
        """make sure the reads to be evaluated actually span the site of the cut."""
        check = True
        for read in reads:
            amplicon = read.reference_name
            interval = self.event_window_by_amplicon[amplicon]
            reference_positions = read.get_reference_positions()
            min_pos = np.min(reference_positions)
            max_pos = np.max(reference_positions)
            check = check and (min_pos < min(interval) & max_pos > max(interval))
        return check

    def extract_events(self, reads: Tuple[AlignedSegment], counts: int):
        """Turn a read or read pair into a set of events."""
        cigars = []
        aln_pos_start = []
        query_seq = []
        real_ref = []
        read_id = []
        mapq = []
        seqnames = []
        strands = []
        for i, r in enumerate(reads, 1):
            reference_name = r.reference_name
            cigars.append(r.cigarstring)
            aln_pos_start.append(r.reference_start + 1)
            real_ref.append(self._ref_lookup[reference_name])
            read_id.append(f"R{i}")
            mapq.append(r.mapping_quality)
            seqnames.append(reference_name)
            if r.is_reverse:
                strands.append("-")
                query_seq.append(
                    mbf_genomes.common.reverse_complement(r.get_forward_sequence())
                )
            else:
                strands.append("+")
                query_seq.append(r.get_forward_sequence())
        events = self.evaluate_cigar(
            cigars,
            aln_pos_start,
            query_seq,
            real_ref,
            read_id,
            mapq,
            seqnames,
            strands,
            [counts] * len(cigars),
        )
        return events

    def map_events_to_reference_space(self, df: DataFrame, reads: List[AlignedSegment]):
        """Convert all event coordinates to reference space"""
        ref_start, ref_end = [], []
        for _, row in df.iterrows():
            read = reads[int(row["read_id"][1])-1]
            start, end = self.map_pairwise_space_to_reference_space(
                read, row["start"], row["end"], row["type"]
            )
            ref_start.append(start)
            ref_end.append(end)
        df["start"] = ref_start
        df["end"] = ref_end
        return df

    def filter_events(
        self, events: DataFrame, is_paired: bool, reads: Tuple[AlignedSegment]
    ):
        """Filter events"""
        keys = [
            "_".join(
                [
                    str(row[x])
                    for x in ["seqnames", "type", "originally", "replacement", "start", "end"]
                ]
            )
            for _, row in events.iterrows()
        ]
        events["keys"] = keys
        filtered = []
        # check pe reads only if both reads agree or a single read is sufficient
        support = 1
        if is_paired and (not self.single_read_sufficient):
            support = 2
        for _, df_sub in events.groupby("keys"):
            # decide what to do
            amplicon = df_sub["seqnames"].values[0]
            interval = self.ranges_by_amplicon[amplicon]
            # throw out events beyond the range of interest
            if (
                df_sub["start"].values[0] < interval[0]
                or df_sub["start"].values[0] > interval[1]
            ):
                continue
            event_type = df_sub["type"].values[0]
            if event_type in self.types_to_consider:
                if len(df_sub) == support:
                    ev = df_sub.iloc[[0]]
                    filtered.append(ev)
                else:
                    pass  # discard as the reads do not agree
        res = pd.concat(filtered)
        return res

    def tally_reads(self, sample: Sample):
        """group all mate pairs together and tally the number of idential reads"""

        def __count():
            tmp_dict = {}
            with pysam.AlignmentFile(sample.bam_filename, "rb") as sam:
                for read in sam.fetch():  # only mapped reads
                    if read.is_unmapped:
                        if read.query_name in tmp_dict:
                            tmp_dict[read.query_name].append(read)
                        else:
                            tmp_dict[read.query_name] = [read]
            counter = collections.Counter()
            read_dict = {}
            for name in tmp_dict:
                rs = tmp_dict[name]
                if self.paired and len(rs) != 2:
                    continue
                seqs = []
                for r in rs:
                    seqs.append(r.get_forward_sequence())
                key = tuple(seqs)
                counter[key] += 1
                read_dict[key] = tuple(rs)  # we need only one pair if we count them
            return counter, read_dict

        return __count()

    def count_events(self, sample: Sample, result_dir: Path = None):
        """
        Creates two tables, first one includes a per-read event count, the second one counts each event separately.
        """
        if result_dir is None:
            result_dir = sample.result_dir
        else:
            result_dir = result_dir
        outfile = result_dir / f"Pycheck_{sample.name}_per_read.tsv"
        outfile2 = result_dir / f"Pycheck_{sample.name}_per_event.tsv"
        outfile3 = result_dir / f"Pycheck_{sample.name}_tallied_reads.tsv"

        def __count():
            per_read = {
                "Sample": [],
                "Reference": [],
                "Type": [],
                "Count": [],
                "Info": [],
            }
            per_event_columns = [
                "Sample",
                "Reference",
                "Start",
                "Stop",
                "Type",
                "Original",
                "Replacement",
            ]
            counter, read_dict = self.tally_reads(sample)
            counter_events = collections.Counter()
            with outfile3.open("w") as outp:
                outp.write("Sequence\tCount\n")
                for seq_tuple in read_dict:
                    reads = read_dict[seq_tuple]
                    for r in reads:
                        if r.cigarstring is None or r.cigarstring == "":
                            continue
                    counts = counter[seq_tuple]
                    outp.write(",".join(list(seq_tuple)) + f"\t{counts}\n")
                    events = self.extract_events(reads, counts)
                    if events is None:
                        continue
                    # remap event coordinates to reference space
                    events = self.map_events_to_reference_space(events, reads)
                    events = self.filter_events(events, sample.is_paired, reads)
                    if len(events) == 0:
                        continue
                    info = []
                    for _, row in events.iterrows():
                        event_info = [sample.name] + [
                            str(row[x])
                            for x in [
                                "seqnames",
                                "start",
                                "end",
                                "type",
                                "originally",
                                "replacement",
                            ]
                        ]
                        counter_events[tuple(event_info)] += 1
                        info.append(
                            f"(type:{row['type']},start={row['start']},end={row['end']},original={row['originally']},replacement={row['replacement']})"
                        )
                    if len(events) > 1:
                        event_type = "multi"
                    else:
                        event_type = events["type"].values[0]
                    per_read["Sample"].append(sample.name)
                    per_read["Reference"].append(events["seqnames"].values[0])
                    per_read["Type"].append(event_type)
                    per_read["Count"].append(events["counts"].values[0])
                    per_read["Info"].append(";".join(info))
                per_read = pd.DataFrame(per_read)
                per_event = pd.DataFrame(
                    list(counter_events.keys()), columns=per_event_columns
                )
                per_event["Counts"] = list(counter_events.values())
                per_read.to_csv(outfile, sep="\t", index=False)
                per_event.to_csv(outfile2, sep="\t", index=False)
        return ppg.MultiFileGeneratingJob([outfile, outfile2, outfile3], __count).depends_on(
            sample.load()
        ).depends_on(self.get_dependencies()).depends_on(self.__set_primer_range())

    def map_pairwise_space_to_reference_space(
        self, read: AlignedSegment, start: int, end: int, eventtype: str
    ):
        """Maps the pairwise space coordinates to reference space for insertions"""
        if eventtype == "insertion":
            alinged_pairs = read.get_aligned_pairs()
            try:
                ref_start = alinged_pairs[start - 1][1]
            except IndexError:
                if start == 0:
                    ref_start = read.reference_start - 1
                else:
                    raise
            try:
                ref_end = alinged_pairs[end + 1][1]
            except IndexError:
                if end >= len(alinged_pairs):
                    ref_end = 1 + np.max(read.get_reference_positions())
                else:
                    raise
        elif eventtype == "deletion":
            ref_start = start
            ref_end = end
        elif eventtype == "mismatch":
            ref_start = read.get_aligned_pairs()[start][1]
            ref_end = ref_start        
        else:
            raise NotImplementedError()
        return ref_start, ref_end

    def evaluate_cigar(
        self,
        extendedcigar: List[str],
        aln_pos_start: List[int],
        query_seq: List[str],
        ref: List[str],
        read_id: List[str],
        mapq: List[int],
        seqnames: List[str],
        strands: List[str],
        counts: List[int],
    ):
        ro.r("library(amplican)")
        events_granges = ro.r(
            """
        cigarsToEvents <- function(cigars, aln_pos_start, query_seq, ref, read_id, mapq,
                           seqnames, strands, counts) {
            if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
                stop("Install GenomicAlignments before calling this function.")
            }

            ids <- seq_along(cigars)
            # INS
            ins <- GenomicAlignments::cigarRangesAlongQuerySpace(cigars, ops = "I")
            repl <- Biostrings::extractAt(DNAStringSet(query_seq), ins)
            ins <- GenomicAlignments::cigarRangesAlongPairwiseSpace(cigars, ops = "I")
            #csum <- lapply(ins, function(x) -1L * cumsumw(x))
            #csum <- IRanges::IntegerList(csum)

            names(ins) <- ids
            #names(csum) <- ids
            #csum <- unlist(csum, use.names = TRUE)
            ins <- unlist(ins, use.names = TRUE) # empty ranges are droped out
            #csum <- csum[names(csum) %in% names(ins)]
            #ins <- IRanges::shift(ins, csum) # shift by cumsum of ins

            iids <- as.integer(names(ins))
            ins <- if (length(ins) > 0) {
                GenomicRanges::GRanges(
                seqnames = seqnames[iids],
                ranges = ins,
                strand = strands[iids],
                originally = "",
                replacement = as.character(unlist(repl, use.names = FALSE)),
                type = "insertion",
                read_id = read_id[iids],
                score = mapq[iids],
                counts = counts[iids])
            } else {
                GenomicRanges::GRanges()
            }
            # DEL
            del <- GenomicAlignments::cigarRangesAlongReferenceSpace(
                cigars, ops = c("D", "N"), pos = aln_pos_start)
            origin <- Biostrings::extractAt(DNAStringSet(ref), del)
            names(del) <- ids
            del <- unlist(del, use.names = TRUE)
            iids <- as.integer(names(del))

            del <- if (length(del) > 0) {
                GenomicRanges::GRanges(
                seqnames = seqnames[iids],
                ranges = del,
                strand = strands[iids],
                originally = as.character(unlist(origin, use.names = FALSE)),
                replacement = "",
                type = "deletion",
                read_id = read_id[iids],
                score = mapq[iids],
                counts = counts[iids])
            } else {
                GenomicRanges::GRanges()
            }

            # MISMATCH - X #TODO MD tags contain mm too if X not present in CIGARS
            mm <- GenomicAlignments::cigarRangesAlongQuerySpace(cigars, ops = "X")
            repl <- Biostrings::extractAt(DNAStringSet(query_seq), mm)
            repl <- unlist(repl, use.names = FALSE)
            repl_l <- unlist(strsplit(as.character(repl[width(repl) > 1]), ""))
            repl <- as.character(repl[width(repl) == 1])
            # mm <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigars, ops = "X")
            mm <- GenomicAlignments::cigarRangesAlongPairwiseSpace(cigars, ops = "X")
            names(mm) <- ids
            mm <- unlist(mm, use.names = TRUE)
            mm_l <- mm[width(mm) > 1]
            mm <- mm[width(mm) == 1]
            mm_ln <- names(mm_l)
            mm_l <- IRanges::tile(mm_l, width = 1L)
            names(mm_l) <- mm_ln
            mm_l <- unlist(mm_l, use.names = TRUE)
            iids <- as.integer(c(names(mm), names(mm_l)))

            mm <- if (length(mm) > 0) {
                GenomicRanges::GRanges(
                seqnames = seqnames[iids],
                ranges = c(mm, mm_l),
                strand = strands[iids],
                originally = "",
                replacement = c(repl, repl_l),
                type = "mismatch",
                read_id = read_id[iids],
                score = mapq[iids],
                counts = counts[iids])
            } else {
                GenomicRanges::GRanges()
            }
            seqnames <- unique(as.character(seqnames))
            GenomeInfoDb::seqlevels(mm) <- seqnames
            GenomeInfoDb::seqlevels(del) <- seqnames
            GenomeInfoDb::seqlevels(ins) <- seqnames
            ret = c(mm, del, ins)
            if (length(ret) == 0){
                ret = NULL
            }
        }
        """
        )(
            ro.vectors.StrVector(extendedcigar),
            ro.vectors.IntVector(aln_pos_start),
            ro.vectors.StrVector(query_seq),
            ro.vectors.StrVector(ref),
            ro.vectors.StrVector(read_id),
            ro.vectors.IntVector(mapq),
            ro.vectors.StrVector(seqnames),
            ro.vectors.StrVector(strands),
            ro.vectors.IntVector(counts),
        )
        if events_granges == ro.rinterface.NULL:
            return None
        events_df = ro.r("function(events){as.data.frame(events, row.names=c(1:length(events)))}")(events_granges)
        events_df = mbf_r.convert_dataframe_from_r(events_df)
        events_df = events_df.astype({"seqnames": "str"})
        events_df = events_df.reset_index()
        # adjust for 0-based genome index
        events_df["start"] = events_df["start"] - 1
        events_df["end"] = events_df["end"] - 1
        return events_df
