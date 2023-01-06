#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""amplican.py: Contains ...."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from pypipegraph import Job
from mbf.align import Sample
from pandas import DataFrame
from collections import defaultdict
import pandas as pd
import pypipegraph as ppg
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as numpy2ri
import time
import pysam

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class Amplican:
    def __init__(
        self,
        name: str = None,
        min_freq: float = 0.01,
        normalize_by=["guideRNA", "Group"],
        **kwargs,
    ) -> None:
        """
        Wrapper class for the AmpliCan tool.

        This uses AmpliCan to analyze genome editing events.

        https://www.bioconductor.org/packages/release/bioc/html/amplican.html

        The tool is initalized with default parameters which can all be supplied
        via kwargs as required.
        The idea is to use this on all demultiplexed lanes in one go.
        In addition, a dataframe must be supplied, that resembles the config file
        needed for AmpliCan. Since this varies strongly for each experiment,
        this is not included here.
        Adjust for potential library bleed via min_freq, default is 0.01.

        Parameters
        ----------
        min_freq : float
            Minimum event frequency, default is 1 percent. Figure
        normalize_by : list, optional
            Name of config columns to normalize by, by default ["guideRNA", "Group"].
            If no control given, no normalization will be performed.
        """
        if name is None:
            self.name = "Amplican"
        else:
            self.name = name
        self.average_quality = kwargs.get("average_quality", 0)
        self.min_quality = kwargs.get("min_quality", 0)
        self.use_parallel = True
        self.gap_opening = kwargs.get("gap_opening", 25)
        self.gap_extension = kwargs.get("gap_extension", 0)
        self.fastqfiles = kwargs.get("fastqfiles", 0.5)
        self.primer_mismatch = kwargs.get("primer_mismatch", 0)
        self.donor_mismatch = kwargs.get("donor_mismatch", 1)
        self.primer_dimer = kwargs.get("PRIMER_DIMER", 30)
        self.event_filter = kwargs.get("event_filter", True)
        self.cut_buffer = kwargs.get("cut_buffer", 5)
        self.promiscuous_consensus = kwargs.get("promiscuous_consensus", True)
        self.min_freq = min_freq
        self.cache_dir = Path("cache") / self.name
        self.normalize = normalize_by
        self.parallel = True

    def get_dependencies(self) -> List[Job]:
        """
        Returns the list of job dependencies.

        Returns
        -------
        List[Job]
            List of job dependencies.
        """
        return [
            ppg.ParameterInvariant(
                f"{self.name}_init_params",
                [
                    self.min_quality,
                    self.average_quality,
                    self.gap_opening,
                    self.fastqfiles,
                    self.primer_mismatch,
                    self.donor_mismatch,
                    self.primer_dimer,
                    self.event_filter,
                    self.cut_buffer,
                    self.promiscuous_consensus,
                    self.min_freq,
                    self.normalize,
                ],
            )
        ]

    def prepare_config(
        self, df_loading_function: Callable, outfile: Path, dependencies: List[Job] = []
    ) -> Job:
        """
        Writes the config file for a run call.

        Parameters
        ----------
        df_loading_function : Callable
            DataFrame loading function.
        outfile : Path
            File path to write to.
        dependencies : List[Job], optional
            List of dependencies, by default [].

        Returns
        -------
        Job
            The job writing the file.
        """
        outfile.parent.mkdir(parents=True, exist_ok=True)
        dependencies.append(ppg.FunctionInvariant(f"prepare_config_{outfile}", df_loading_function))

        def __dump():
            df = df_loading_function()
            df.to_csv(outfile, index=False)

        return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)

    def prepare_fastqs(self, samples: Dict[str, Sample], fastq_folder: Path) -> List[Job]:
        """
        Returns a list of dependencies for the run job.

        Since we need all input files in the same directory, this modifies the
        prepare_input function of the samples in order to produces named
        files in a single directory. This is achieved by modifying the
        get_aligner_input_files function.

        Parameters
        ----------
        samples : Dict[str, Sample]
            Dictionary of raw samples.
        fastq_folder : Path
            The common output folder.

        Returns
        -------
        List[Job]
            List of job dependencies self.run needs to depend on.
        """

        def common_path_filenames(sample):
            def __get_aligner_input_filenames():
                if sample.is_paired:
                    return (
                        fastq_folder / f"{sample.name}_R1_.fastq",
                        fastq_folder / f"{sample.name}_R2_.fastq",
                    )
                else:
                    return (fastq_folder / f"{sample.name}_R1_.fastq",)

            return __get_aligner_input_filenames

        dependencies = []
        for sample_name in samples:
            sample = samples[sample_name]
            sample.get_aligner_input_filenames = common_path_filenames(sample)
            dependencies.append(sample.prepare_input())
        return dependencies

    def run(
        self,
        run_id: str,
        samples: Dict[str, Sample],
        df_config_job: Job,
        output_path: Union[str, Path],
        dependencies: List[Job] = [],
    ) -> Job:
        """
        Runs AmpliCan for all lanes given.

        Since the amplican run needs a config file for each run, a dataframe
        containing the following columns must be supplied:
            ID, Barcode, Forward_Reads, Reverse_Reads, Group, Control, guideRNA,
            Forward_Primer, Reverse_Primer, Direction, Amplicon, Donor.
        This is handled by the df_config_job, since this will change with each experiment.

        Parameters
        ----------
        run_id : str
            A unique name for the run.
        samples : Dict[str, Sample]
            Raw samples to analyze.
        df_config_job : Job
            Job that writes the config DataFrame.
        output_path : Union[str, Path]
            Result directory.
        dependencies : List[Job], optional
            Additional dependencies, by default [].

        Returns
        -------
        Job
            Job that runs the analysis and creates a sentinel file.
        """
        fastq_folder = self.cache_dir / run_id / "fastq"
        fastq_folder.mkdir(exist_ok=True, parents=True)
        input_dependencies = self.prepare_fastqs(samples, fastq_folder)
        if isinstance(output_path, str):
            result_dir = Path(output_path)
        else:
            result_dir = output_path
        result_dir.mkdir(parents=True, exist_ok=True)
        outfile = result_dir / "sentinel.txt"

        def __dump():
            with outfile.open("w") as op:
                start_time = time.time()
                self.call_amplican(result_dir, fastq_folder, str(df_config_job.job_id))
                runtime = time.time() - start_time
                op.write(f"run time: {runtime:.2f} seconds")

        job = ppg.FileGeneratingJob(outfile, __dump)
        job.depends_on(df_config_job)
        job.depends_on(input_dependencies)
        job.depends_on(dependencies)
        return job

    def call_amplican(self, result_dir: Path, fastq_folder: Path, config_file: str):
        """
        Actually calls the amplican method via r2py.

        Parameters
        ----------
        result_dir : Path
            The oputput directory.
        fastq_folder : Path
            the common fastq folder.
        config_file : str
            The config file.
        """
        print(config_file)
        print(str(fastq_folder))
        print(str(result_dir))
        ro.r("library(amplican)")
        ro.r("amplicanPipeline")(
            config=config_file,
            fastq_folder=str(fastq_folder),
            results_folder=str(result_dir),
            knit_reports=True,
            write_alignments_format="txt",
            average_quality=self.average_quality,
            min_quality=self.min_quality,
            use_parallel=self.parallel,
            # scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(match = 5,mismatch = -4, baseOnly = TRUE, type = "DNA"), use default for now
            gap_opening=self.gap_opening,
            gap_extension=self.gap_extension,
            fastqfiles=self.fastqfiles,
            primer_mismatch=self.primer_mismatch,
            donor_mismatch=self.donor_mismatch,
            PRIMER_DIMER=self.primer_dimer,
            event_filter=self.event_filter,
            cut_buffer=self.cut_buffer,
            promiscuous_consensus=self.promiscuous_consensus,
            normalize=ro.vectors.StrVector(self.normalize),
            min_freq=self.min_freq,
        )

    @classmethod
    def cigar_to_event(
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
        if counts is None:
            counts = ro.NULL
        events_granges = ro.r("cigarsToEvents")(
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
        events_df = ro.r("function(events){as.data.frame(events, row.names=c(1:length(events)))}")(
            events_granges
        )
        events_df = mbf_r.convert_dataframe_from_r(events_df)
        events_df = events_df.astype({"seqnames": "str"})
        events_df = events_df.reset_index()
        return events_df

    @classmethod
    def write_top_x_alignments(
        self,
        output_file: Union[str, Path],
        result_dir: Path,
        x: int = 50,
        dependencies: List[Job] = [],
    ):
        """
        Gets the x most common alignments per sample and writes it to file.

        Parameters
        ----------
        output_file : Union[str, Path]
            Output file name or path.
        result_dir : Path
            Amplican result dir, where the alignment folder is located.
        x : int, optional
            Number of most common alignments to collect, by default 50.
        dependencies : List[Job]
            Job dependencies.
        """
        if isinstance(output_file, str):
            outfile = Path(output_file)
        else:
            outfile = output_file
        outfile.parent.mkdir(parents=True, exist_ok=True)

        def __write():
            return self.get_top_x_alignments(result_dir, x)

        return ppg.FileGeneratingJob(outfile, __write).depends_on(dependencies)

    @classmethod
    def get_top_x_alignments(self, result_dir: Path, x: int):
        """
        Gets the x most common alignments per sample and returns it in a
        dataframe.

        Parameters
        ----------
        result_dir : Path
            Amplican result dir, where the alignment folder is located.
        x : int, optional
            Number of most common alignments to collect, by default 50.
        """
        to_df: Dict[str, List] = {
            "ID": [],
            "Count": [],
            "Sequence": [],
            "Reference": [],
        }
        current_id = None
        with (result_dir / "alignments" / "alignments.txt").open("r") as inp:
            for next_id, count, seq, ref in next_4(inp):
                if next_id != current_id:
                    counter = 0
                    current_id = next_id
                if counter < x:
                    to_df["ID"].append(next_id)
                    to_df["Count"].append(count)
                    to_df["Sequence"].append(seq)
                    to_df["Reference"].append(ref)
        return pd.DataFrame(to_df)


"""
def create_bam_from_alignments(sample_name, output_file: Union[str, Path], result_dirs: Dict[str, Path], genome: Genome, dependencies: List[Job] = []):
    if isinstance(output_file, str):
        outfile = Path(output_file)
    else:
        outfile = output_file
    outfile.parent.mkdir(parents=True, exist_ok=True)

    def get_aligned_segment(reference_name: str, reference_id: str, chrm_name: str, next_id, count, seq, ref):
        segment = pysam.AlignedSegment()
        segment.query_name = f"Alignment_{next_id}_{no}_{count}"
        segment.query_sequence = seq.replace("-", "")
        segment.reference_id = 0
        
        a.flag = 99
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((0,10), (2,1), (0,25))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
    a.tags = (("NM", 1),
              ("RG", "L1"))
    def build_pysam():
        # make header
        sq = []
        for chrname, len in genome.get_chromosome_lengths():
            sq.append({"LN": len, "SN": chrname})
        header = {'HD': {'VN': '1.0'}, 'SQ': sq}
        for chrmname in result_dirs:
            result_dir = result_dirs[chrmname]
            with (result_dir / "alignments" / "alignments.txt").open("r") as inp:
        
        
        
        with pysam.AlignmentFile(outfile, "wb", header=header) as outf:
                for no, (next_id, count, seq, ref) in enumerate(next_4(inp)):
                    if next_id != sample_name:
    outf.write(a)

    return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)    
"""


def next_4(inp):
    row1 = inp.readline()
    row2 = inp.readline()
    row3 = inp.readline()
    inp.readline()
    while row1:
        splits = row1[-1].split()
        print(splits)
        next_id = splits[1]
        count = splits[5]
        ref = row2[:-1]
        seq = row3[:-1]
        yield (next_id, count, seq, ref)
        row1 = inp.readline()
        row2 = inp.readline()
        row3 = inp.readline()
        inp.readline()
