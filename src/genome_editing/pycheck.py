#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pysam.py: Contains ...."""

from typing import List, Dict, Tuple, Union
from mbf.align import Sample

# from mbf_genomes import GenomeBase
from pysam import AlignedSegment
from pandas import DataFrame
from pypipegraph import Job
from pathlib import Path
import pandas as pd
import pypipegraph as ppg
import pysam
import numpy as np

# import mbf_genomes
# import mbf_genomics
# import collections
# import mbf_r
import mbf
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as numpy2ri
import subprocess
import tempfile
import os

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class PysamCheck:
    """The idea is to use pysam to check the cigar strings of alinged reads to
    extract the editing event"""

    # TODO: get the modified sequence --> by drawing the pairwise space numbers
    # TODO: ignore all events outside the range

    def __init__(
        self,
        primers_or_range_by_amplicon_id: Dict[str, List[str]],
        event_window_by_amplicon: Dict[str, List[int]],
        genome: GenomeBase,
        paired: bool = True,
        mapping_qality_filter: int = 20,
        single_read_sufficient: bool = False,
        types_to_consider: List[str] = None,
        name: str = None,
    ):
        self.name = "PysamCheck"
        if name is not None:
            self.name = name
        self.primers_or_range_by_amplicon = (
            primers_or_range_by_amplicon_id  # events within primer ranges are not reported
        )
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
            ppg.FunctionInvariant(
                "map_pairwise_space_to_reference_space",
                self.map_pairwise_space_to_reference_space,
            ),
            ppg.FunctionInvariant("evaluate_cigar", self.evaluate_cigar),
            ppg.FunctionInvariant("tally_reads", self.tally_reads),
            ppg.FunctionInvariant(
                "map_events_to_reference_space", self.map_events_to_reference_space
            ),
            ppg.FunctionInvariant("filter_events", self.filter_events),
            ppg.FunctionInvariant("extract_events", self.extract_events),
            ppg.FunctionInvariant(
                "__check_reads_spanning_windows", self.__check_reads_spanning_windows
            ),
            ppg.FunctionInvariant("__set_primer_range", self.__set_primer_range),
            ppg.FunctionInvariant("_prepare_reference_lookup", self._prepare_reference_lookup),
            self.genome.download_genome(),
        ]
        return deps

    def _prepare_reference_lookup(self):
        """put the reference sequences in a dict for quick lookup"""

        def __do_lookup():
            lookup: Dict[str, str] = {}
            for name, length in self.genome.get_chromosome_lengths().items():
                lookup[name] = self.genome.get_genome_sequence(name, 0, length)
            return lookup

        job = ppg.AttributeLoadingJob(
            f"{self.name}__lookup", self, "_ref_lookup", __do_lookup
        ).depends_on(self.get_dependencies())
        return job

    def __set_primer_range(self):
        """set the range in which we report events"""

        def __set_range_from_primers():
            if isinstance(next(iter(self.primers_or_range_by_amplicon.values()))[0], str):
                ranges: Dict[str, Tuple(int)] = {}
                for amplicon in self.primers_or_range_by_amplicon:
                    primer_forward, primer_reverse = self.primers_or_range_by_amplicon[amplicon]
                    seq = self._ref_lookup[amplicon]
                    i1 = seq.find(primer_forward) + len(primer_forward)
                    i2 = seq.find(primer_reverse) + len(primer_reverse)
                    ranges[amplicon] = (i1, i2)
                return ranges
            else:
                return self.primers_or_range_by_amplicon

        job = ppg.AttributeLoadingJob(
            f"{self.name}__ranges_by_amplicon",
            self,
            "ranges_by_amplicon",
            __set_range_from_primers,
        )
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

    def extract_events(self, reads: Tuple[AlignedSegment], counts: int) -> Union[DataFrame, None]:
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
                query_seq.append(mbf_genomes.common.reverse_complement(r.get_forward_sequence()))
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
            [counts],  #  * len(cigars),
        )
        return events

    def map_events_to_reference_space(self, df: DataFrame, reads: List[AlignedSegment]):
        """Convert all event coordinates to reference space"""
        ref_start, ref_end = [], []
        query_start = []
        query_end = []
        for _, row in df.iterrows():
            read = reads[int(row["read_id"][1]) - 1]
            start, end = self.map_pairwise_space_to_reference_space(
                read, row["start"], row["end"], row["type"]
            )
            if row["type"] == "deletion":
                qstart = -1
                qend = -1
            else:
                qstart, qend = self.map_pairwise_space_to_query_space(
                    read, row["start"], row["end"], row["type"]
                )
            query_start.append(qstart)
            query_end.append(qend)
            ref_start.append(start)
            ref_end.append(end)
        df["start"] = ref_start
        df["end"] = ref_end
        df["qstart"] = query_start
        df["qend"] = query_end
        return df

    def filter_events(self, events: DataFrame, is_paired: bool, reads: Tuple[AlignedSegment]):
        """Filter events"""
        keys = [
            "_".join(
                [
                    str(row[x])
                    for x in [
                        "seqnames",
                        "type",
                        "originally",
                        "replacement",
                        "start",
                        "end",
                    ]
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
            if df_sub["start"].values[0] < interval[0] or df_sub["start"].values[0] > interval[1]:
                continue
            event_type = df_sub["type"].values[0]
            if event_type in self.types_to_consider:
                if len(df_sub) == support:
                    ev = df_sub.iloc[[0]]
                    filtered.append(ev)
                else:
                    pass  # discard as the reads do not agree
        if len(filtered) > 0:
            res = pd.concat(filtered)
            return res
        else:
            return None

    def tally_reads(self, sample: Sample):
        """
        Group all mate pairs together and tally the number of idential reads.
        Orphan reads or such with unusable cigar string are discarded.
        """

        def __count():
            tmp_dict = {}
            with pysam.AlignmentFile(sample.bam_filename, "rb") as sam:
                for read in sam.fetch():  # only mapped reads
                    if read.query_name in tmp_dict:
                        tmp_dict[read.query_name].append(read)
                    else:
                        tmp_dict[read.query_name] = [read]
            counter = collections.Counter()
            read_dict = {}
            for name in tmp_dict:
                rs = tmp_dict[name]
                if self.paired and len(rs) != 2:
                    continue  # no orphans
                not_usable = any([(r.cigarstring is None or r.cigarstring == "") for r in rs])
                if not_usable:
                    continue
                key = tuple([r.get_forward_sequence() for r in rs])
                counter[key] += 1
                read_dict[key] = tuple(rs)  # we need only one pair if we count them
            return counter, read_dict

        return __count()

    def tally(self, sample: Sample, result_dir: Path = None):
        if result_dir is None:
            result_dir = sample.result_dir
        else:
            result_dir = result_dir
        result_dir.mkdir(exist_ok=True, parents=True)
        outfile3 = result_dir / f"Pycheck_{sample.name}_tallied_reads.tsv"
        outfile4 = result_dir / f"Pycheck_{sample.name}_tallied_reads.bam"
        outfile5 = result_dir / f"Pycheck_{sample.name}_tallied_reads_ns.bam"

        def __count():
            counter, read_dict = self.tally_reads(sample)
            print(sample.bam_filename)
            samfile = pysam.AlignmentFile(sample.bam_filename, "rb")
            print(samfile.header)
            running = 0
            with tempfile.NamedTemporaryFile() as tmp:
                tallied_bam = pysam.AlignmentFile(tmp, "wb", header=samfile.header)
                with outfile3.open("w") as outp:
                    outp.write("Sequence\tCount\n")
                    for order, (seq_tuple, counts) in enumerate(counter.most_common()):
                        reads = read_dict[seq_tuple]
                        #                        counts = counter[seq_tuple]
                        running += counts
                        for read in reads:
                            # add count as a tag
                            read.set_tag("RC", counts, "i")
                            read.set_tag("OD", order, "i")
                            # write the reads
                            tallied_bam.write(read)
                        outp.write(",".join(list(seq_tuple)) + f"\t{counts}\n")
                tallied_bam.close()
                pysam.sort("-o", str(outfile4), tmp.name)
                pysam.index(str(outfile4))
                pysam.sort("-n", "-o", str(outfile5), tmp.name)
            print(running)

        #            raise ValueError()
        return (
            ppg.MultiFileGeneratingJob([outfile3, outfile4], __count)
            .depends_on(sample.load())
            .depends_on(self.get_dependencies())
            .depends_on(self.__set_primer_range())
        )

    def _get_homology(self, deletion: str, flanking_forward: str, flanking_reverse: str):
        hom1 = os.path.commonprefix([deletion, flanking_forward])
        hom2 = os.path.commonprefix([deletion[::-1], flanking_reverse[::-1]])[::-1]
        homologies = [hom1, hom2]
        i = np.argmax(homologies)
        longest_homology = homologies[i]
        return longest_homology

    def write_microhomologies(self, sample: Sample, result_dir: Path):
        result_dir.mkdir(exist_ok=True, parents=True)
        outfile = result_dir / f"Pycheck_{sample.name}_microhomologies.tsv"
        infile = result_dir / f"Pycheck_{sample.name}_tallied_reads_ns.bam"
        infile_per_read = result_dir / f"Pycheck_{sample.name}_per_read.tsv"

        def __write():
            to_df = {
                "Read id": [],
                "Deletion": [],
                "Microhomology": [],
                "Count": [],
                "Reference": [],
            }
            read_dict = {}
            with pysam.AlignmentFile(infile, "rb") as tallied:
                iterator = tallied.fetch(until_eof=True)
                for read in iterator:
                    reads = read
                    if sample.is_paired:
                        mate = iterator.__next__()  # tallied.fetch(until_eof=True)
                        assert read.get_tag("OD") == mate.get_tag("OD")
                        reads = (read, mate)
                        read_dict[read.query_name] = reads
            df_allread = pd.read_csv(infile_per_read, sep="\t")
            for reference, df_read in df_allread.groupby("Reference"):
                for _, row in df_read.iterrows():
                    event_type = row["Type"]
                    if event_type != "deletion":
                        continue
                    events = row["Info"].split(",")
                    if len(events) > 1:
                        raise ValueError("Deletion should only contain 1 depetion event.")
                    _, del_start, del_end, original, _, _, _ = events[0].split("_")
                    read_id = row["Read id"]
                    reads = read_dict[read_id]
                    del_start = int(del_start)
                    del_end = int(del_end) + 1
                    after_deletion = self._ref_lookup[row["Reference"]][del_end:]
                    before_deletion = self._ref_lookup[row["Reference"]][:del_start]
                    homology = self._get_homology(original, after_deletion, before_deletion)
                    counts = reads[0].get_tag("RC")
                    to_df["Read id"].append(read_id)
                    to_df["Deletion"].append(original)
                    to_df["Microhomology"].append(homology)
                    to_df["Count"].append(counts)
                    to_df["Reference"].append(reference)
            df = pd.DataFrame(to_df)
            df.to_csv(outfile, sep="\t", index=False)

        dep = [
            self.tally(sample, result_dir),
            self.count_events(sample, result_dir),
            self.__set_primer_range(),
        ]
        return ppg.FileGeneratingJob(outfile, __write).depends_on(dep)

    def count_events(self, sample: Sample, result_dir: Path = None):
        """
        Creates two tables, first one includes a per-read event count, the second one counts each event separately.
        """
        if result_dir is None:
            result_dir = sample.result_dir
        else:
            result_dir = result_dir
        result_dir.mkdir(exist_ok=True, parents=True)
        outfile = result_dir / f"Pycheck_{sample.name}_per_read.tsv"
        outfile2 = result_dir / f"Pycheck_{sample.name}_per_event.tsv"
        infile5 = result_dir / f"Pycheck_{sample.name}_tallied_reads_ns.bam"

        def __count():
            per_read = {
                "Sample": [],
                "Read id": [],
                "Reference": [],
                "Type": [],
                "Count": [],
                "Info": [],
                "Read sequence": [],
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
            counter_events = collections.Counter()
            no_events = 0
            filtered_out = 0
            with pysam.AlignmentFile(infile5, "rb") as tallied:
                iterator = tallied.fetch(until_eof=True)
                for read in iterator:
                    reads = read
                    if sample.is_paired:
                        mate = iterator.__next__()  # tallied.fetch(until_eof=True)
                        assert read.get_tag("OD") == mate.get_tag("OD")
                        reads = (read, mate)
                    counts = read.get_tag("RC")
                    events = self.extract_events(reads, counts)
                    if events is None:
                        no_events += counts
                        continue
                    # remap event coordinates to reference space
                    events = self.map_events_to_reference_space(events, reads)
                    filtered = self.filter_events(events, sample.is_paired, reads)
                    if filtered is None:
                        filtered_out += counts
                        continue
                    info = []
                    for _, row in filtered.iterrows():
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
                        counter_events[tuple(event_info)] += counts
                        info.append(
                            "_".join(
                                [
                                    row["type"],
                                    str(row["start"]),
                                    str(row["end"]),
                                    row["originally"],
                                    row["replacement"],
                                    str(row["qstart"]),
                                    str(row["qend"]),
                                ]
                            )
                        )
                    if len(filtered) > 1:
                        event_type = "multi"
                    else:
                        event_type = filtered["type"].values[0]
                    read_seq = []
                    for read_id, filtered_part in filtered.groupby("read_id"):
                        qstarts = [x for x in filtered_part["qstart"].values if x >= 0]
                        r_start = 0 if len(qstarts) == 0 else min(qstarts)
                        qends = [x for x in filtered_part["qend"].values if x >= 0]
                        r_end = 0 if len(qends) == 0 else max(qends) + 1
                        seq_index = int(read_id[1]) - 1
                        read_seq.append(reads[seq_index].get_forward_sequence()[r_start:r_end])
                    per_read["Sample"].append(sample.name)
                    per_read["Reference"].append(filtered["seqnames"].values[0])
                    per_read["Read id"].append(reads[0].query_name)
                    per_read["Type"].append(event_type)
                    per_read["Count"].append(counts)
                    per_read["Info"].append(",".join(info))
                    per_read["Read sequence"].append(",".join(read_seq))
            per_read = pd.DataFrame(per_read)
            per_event = pd.DataFrame(list(counter_events.keys()), columns=per_event_columns)
            per_event["Counts"] = list(counter_events.values())
            # add no events
            row_df = pd.DataFrame(
                [
                    [sample.name, "", "", "", "No event", "", "", no_events],
                    [sample.name, "", "", "", "Filtered out", "", "", filtered_out],
                ],
                columns=per_event_columns + ["Counts"],
            )
            per_event = pd.concat([row_df, per_event], ignore_index=True)
            per_read.to_csv(outfile, sep="\t", index=False)
            per_event.to_csv(outfile2, sep="\t", index=False)

        return (
            ppg.MultiFileGeneratingJob([outfile, outfile2], __count)
            .depends_on(sample.load())
            .depends_on(self.get_dependencies())
            .depends_on(self.__set_primer_range())
            .depends_on(self.tally(sample, result_dir))
        )

    def map_pairwise_space_to_reference_space(
        self, read: AlignedSegment, start: int, end: int, eventtype: str
    ):
        """Maps the pairwise space coordinates to reference space for insertions"""
        offset = 0
        if read.cigartuples[0][0] in [4]:
            offset = read.cigartuples[0][1]
        if eventtype == "insertion":
            alinged_pairs = read.get_aligned_pairs()
            try:
                ref_start = alinged_pairs[start - 1 + offset][1]
            except IndexError:
                if start == 0:
                    ref_start = read.reference_start - 1
                else:
                    raise
            try:
                ref_end = alinged_pairs[end + 1 + offset][1]
            except IndexError:
                if end >= len(alinged_pairs):
                    ref_end = 1 + np.max(read.get_reference_positions())
                else:
                    raise
        elif eventtype == "deletion":
            ref_start = start
            ref_end = end
        elif eventtype == "mismatch":
            ref_start = read.get_aligned_pairs()[start + offset][1]
            ref_end = ref_start
        else:
            raise NotImplementedError()
        return ref_start, ref_end

    def map_pairwise_space_to_query_space(
        self, read: AlignedSegment, start: int, end: int, eventtype: str
    ):
        """Maps the pairwise space coordinates to query space for insertions/mismatches"""
        offset = 0
        if read.cigartuples[0][0] in [4]:
            offset = read.cigartuples[0][1]
        if eventtype == "insertion":
            alinged_pairs = read.get_aligned_pairs()
            ref_start = alinged_pairs[start + offset][0]
            ref_end = alinged_pairs[end + offset][0]
        elif eventtype == "deletion":
            raise ValueError("Deletions are not in the query")
        elif eventtype == "mismatch":
            ref_start = read.get_aligned_pairs()[start + offset][0]
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
    ) -> Union[DataFrame, None]:
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
        events_df = ro.r(
            """
            function(events){
                l = length(events)
                if (l == 0) {
                    NULL
                }
                else {
                    as.data.frame(events, row.names=c(1:l))
                }
            }
            """
        )(events_granges)
        if events_df == ro.r("NULL"):
            return None
        events_df = mbf_r.convert_dataframe_from_r(events_df)
        events_df = events_df.astype({"seqnames": "str"})
        events_df = events_df.reset_index()
        # adjust for 0-based genome index
        events_df["start"] = events_df["start"] - 1
        events_df["end"] = events_df["end"] - 1
        return events_df
