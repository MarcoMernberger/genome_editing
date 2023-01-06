#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import time
import random
import docker
import pypipegraph as ppg
import mbf
import re
import subprocess
from mbf.genomes.ensembl import _EnsemblGenome
from mbf.align.raw import Sample
from pathlib import Path
from typing import List, Any, Union, Optional, Dict
from pandas import DataFrame
from pypipegraph2 import Job
from mutility import reverse_complement, dict_to_string_of_items, read_excel_from_biologists


def create_crispresso_input_frame(
    path_to_input_table: Path,
    genome: _EnsemblGenome,
    amplicon_name_column: str = "Gen",
    primer_fwd_column: str = "Primer forward",
    primer_rev_column: str = "Primer reverse",
    sgRNA_column: str = "shRNA",
    amplicon_column: str = "Amplicon",
    dependencies: List[Job] = [],
):
    outfile = Path("cache/crispresso2/df_amplicon.tsv")
    outfile.parent.mkdir(parents=True, exist_ok=True)
    df_in = read_excel_from_biologists(path_to_input_table)
    rename = {
        amplicon_name_column: "Gen",
        primer_fwd_column: "Primer forward",
        primer_rev_column: "Primer reverse",
        sgRNA_column: "shRNA",
        amplicon_column: "Amplicon",
    }
    df_out = df_in.rename(columns=rename)
    if "Amplicon" not in df_in.columns:
        job = write_amplicon_sequences_by_name(genome, df_out, outfile, dependencies)
    else:

        def __write(outfile):
            df_out.to_csv(outfile, sep="\t", index=False)

        job = ppg.FileGeneratingJob(outfile, __write).depends_on(dependencies)
    return job


def write_amplicon_sequences_by_name(
    genome: _EnsemblGenome,
    df_amplicons: DataFrame,
    output_file: Union[str, Path],
    dependencies: List[Job] = [],
) -> Job:
    outfile = Path(output_file)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    def write_amplicon_df(outfile):
        amplicons = []
        for _, row in df_amplicons.iterrows():
            gene_name = row["Gen"].upper()
            primer_fwd = row["Primer forward"]
            primer_rev = reverse_complement(row["Primer reverse"])
            sgrna = row["sgRNA"]
            if not gene_name.startswith("CHR"):
                ensemblid = genome.name_to_canonical_id(gene_name.upper())
                gene = genome.genes[ensemblid]
                s = genome.get_genome_sequence(gene.chr, gene.start, gene.stop)
                try:
                    amplicon = get_amplicon_sequence(primer_fwd, primer_rev, s, sgrna)
                except ValueError:
                    try:
                        s = get_chromosome_sequence(gene.chr, genome)
                        amplicon = get_amplicon_sequence(primer_fwd, primer_rev, s, sgrna)
                    except ValueError as err:
                        raise ValueError(
                            str(err) + f"\nSearched chromosome {gene.chr} for {gene_name}."
                        )
            else:
                # this is not a gene, but a chromosome
                chr = gene_name.replace("CHR", "")
                s = get_chromosome_sequence(chr, genome)
                amplicon = get_amplicon_sequence(primer_fwd, primer_rev, s, sgrna)
            #            if sgrna not in amplicon:
            #                amplicon = reverse_complement(amplicon)
            #            if sgrna not in amplicon:
            #                get_amplicon_sequences(primer_fwd, primer_rev, s, sgrna)
            #                raise ValueError(f"{gene_name} sgrna {sgrna} not in amplicon:\n {amplicon}")
            amplicons.append(amplicon)
        df_amplicons["Amplicon"] = amplicons
        for _, row in df_amplicons.iterrows():
            try:
                assert row["sgRNA"] in row["Amplicon"]
            except AssertionError:
                print(row["Gen"])
                print("sgRNA: ", row["sgRNA"])
                print("Amplicon: ", row["Amplicon"])
                print("sgRNA in Amplicon?", row["sgRNA"] in row["Amplicon"])
                print("sgRNA rev in Amplicon?", reverse_complement(row["sgRNA"]) in row["Amplicon"])
                raise

        df_amplicons.to_csv(outfile, sep="\t", index=False)

    return (
        ppg.FileGeneratingJob(outfile, write_amplicon_df)
        .depends_on(genome.download_genome())
        .depends_on(dependencies)
    )


def get_amplicon_sequence(primer_fwd: str, primer_rev: str, sequence: str, sgrna: str) -> str:
    primer_fwd_rev = reverse_complement(primer_fwd)
    primer_rev_rev = reverse_complement(primer_rev)
    if sgrna not in sequence:
        sequence = reverse_complement(sequence)
    if sgrna not in sequence:
        raise ValueError(f"sgRNA {sgrna} not in sequence.")

    patterns = [
        f".*{primer_fwd}(?P<amplicon>[ATCG]+){primer_rev}.*",
        f".*{primer_fwd_rev}(?P<amplicon>[ATCG]+){primer_rev}.*",
        f".*{primer_fwd}(?P<amplicon>[ATCG]+){primer_rev_rev}.*",
        f".*{primer_fwd_rev}(?P<amplicon>[ATCG]+){primer_rev_rev}.*",
        f".*{primer_rev}(?P<amplicon>[ATCG]+){primer_fwd}.*",
        f".*{primer_rev_rev}(?P<amplicon>[ATCG]+){primer_fwd}.*",
        f".*{primer_rev}(?P<amplicon>[ATCG]+){primer_fwd_rev}.*",
        f".*{primer_rev_rev}(?P<amplicon>[ATCG]+){primer_fwd_rev}.*",
    ]
    for pattern in patterns:
        match = re.match(pattern, sequence)
        if match is not None:
            break
    if match is None:
        if primer_fwd not in sequence:
            print("primer_fwd in sequence", primer_fwd in sequence)
            print("primer_fwd rev in sequence", reverse_complement(primer_fwd) in sequence)
            print("primer_rev in sequence", primer_rev in sequence)
            print("primer_rev rev in sequence", reverse_complement(primer_rev) in sequence)
            raise ValueError("Primer forward not in sequence.")
        elif primer_rev not in sequence:
            raise ValueError("Primer reverse not in sequence.")
        raise ValueError("No matching patterns found. Are the primers in sequence?")
    else:
        if len(match.groups()) == 1:
            return match.group("amplicon")
        else:
            raise ValueError(f"Multiple possible amplicons detected!\n{match.groups()}.")


def get_amplicon_sequences(primer_fwd: str, primer_rev: str, sequence: str, sgrna: str) -> str:
    primer_fwd_rev = reverse_complement(primer_fwd)
    primer_rev_rev = reverse_complement(primer_rev)
    if sgrna not in sequence:
        sequence = reverse_complement(sequence)
    if sgrna not in sequence:
        raise ValueError(f"This sgrna {sgrna} not in sequence:\n {sequence}")
    patterns = [
        f".*{primer_fwd}(?P<amplicon>[ATCG]+){primer_rev}.*",
        f".*{primer_fwd_rev}(?P<amplicon>[ATCG]+){primer_rev}.*",
        f".*{primer_fwd}(?P<amplicon>[ATCG]+){primer_rev_rev}.*",
        f".*{primer_fwd_rev}(?P<amplicon>[ATCG]+){primer_rev_rev}.*",
        f".*{primer_rev}(?P<amplicon>[ATCG]+){primer_fwd}.*",
        f".*{primer_rev_rev}(?P<amplicon>[ATCG]+){primer_fwd}.*",
        f".*{primer_rev}(?P<amplicon>[ATCG]+){primer_fwd_rev}.*",
        f".*{primer_rev_rev}(?P<amplicon>[ATCG]+){primer_fwd_rev}.*",
    ]
    amplicons = []
    for pattern in patterns:
        matches = re.findall(pattern, sequence)
        amplicons.extend(matches)
    if len(amplicons) == 0:
        raise ValueError("No matching patterns found. Are the primers in sequence?")
    return amplicons


def get_chromosome_sequence(chr, genome):
    lengths = genome.get_chromosome_lengths()
    s = genome.get_genome_sequence(chr, 1, lengths[chr])
    return s


# class Crispresso2:
#     """
#      Wrapper for the crispresso2 docker.

#     [extended_summary]

#     Returns
#     -------
#     [type]
#         [description]
#     """

#     def __init__(self, fastq_dir="/project/incoming"):
#         self.image = "pinellolab/crispresso2"
#         self.wdir = "/project"
#         self.volumes = {
#             os.environ["PWD"]: {"bind": "/project", "mode": "rw"},
#         }
#         self.client = docker.from_env()

#     def get_version(self):
#         s = self.get_help()
#         version = s[s.find("CRISPResso version ") + 19 : s.find("]")]
#         return version

#     def get_help(self):
#         command_help = ["CRISPResso", "-h"]
#         container = self.client.containers.run(
#             image=self.image,
#             volumes=self.volumes,
#             working_dir=self.wdir,
#             command=command_help,
#             detach=True,
#         )
#         ret = ""
#         for line in container.logs(stream=True):
#             ret += line.decode("utf-8")

#         return ret

#     def print_help(self):
#         print(self.get_help())

#     def crispresso_run(
#         self,
#         raw_sample: mbf.align.raw.Sample,
#         genome,
#         report_name: str,
#         df_amplicons: DataFrame,
#         output_folder: Path,
#         dependencies: List[Job] = [],
#         options: List[str] = [],
#         quantification_window_size=10,
#         quantification_window_center=-3,
#     ):
#         if "Amplicon" in df_amplicons:
#             amplicon_seq = ",".join(df_amplicons.Amplicon.values)
#         else:
#             raise NotImplementedError
#         sgrnas = ",".join(df_amplicons.sgRNA.values)
#         amplicon_names = ",".join(df_amplicons.Gen.values)

#         inputfiles = raw_sample.get_aligner_input_filenames()
#         if output_folder is None:
#             output_folder = f"results/crispresso_{self.get_version()}/{raw_sample.name}"
#         output_folder = Path(output_folder)
#         outfile = output_folder / (report_name + ".log")
#         outfile.parent.mkdir(parents=True, exist_ok=True)
#         deps = dependencies
#         deps.append(
#             ppg.ParameterInvariant(
#                 f"PI_{outfile}",
#                 [
#                     self.image,
#                     self.volumes,
#                     self.wdir,
#                     str(options),
#                     quantification_window_center,
#                     quantification_window_size,
#                     amplicon_seq,
#                     sgrnas,
#                     amplicon_names,
#                 ],
#             )
#         )
#         deps.append(raw_sample.prepare_input())

#         def __dump(output_file):
#             command = [
#                 "CRISPResso",
#                 "--fastq_r1",
#                 str(inputfiles[0]),
#                 "--amplicon_seq",
#                 amplicon_seq,
#                 "--guide_seq",
#                 sgrnas,
#                 "--quantification_window_size",
#                 f"{quantification_window_size}",
#                 "--quantification_window_center",
#                 f"{quantification_window_center}",
#                 "--base_editor_output",
#                 "--exclude_bp_from_right",
#                 "1",
#                 "--exclude_bp_from_left",
#                 "1",
#                 "-o",
#                 str(output_folder),
#                 "-n",
#                 report_name,
#                 "-an",
#                 amplicon_names,
#             ]
#             if raw_sample.is_paired:
#                 command.extend(
#                     [
#                         "--fastq_r2",
#                         str(inputfiles[1]),
#                     ]
#                 )
#             for k in options:
#                 if len(str(options[k])) == 0:
#                     command.append(k)
#                 else:
#                     command.extend([k, str(options[k])])
#             time.sleep(random.random() * 2)
#             container = self.client.containers.run(
#                 self.image,
#                 volumes=self.volumes,
#                 working_dir=self.wdir,
#                 command=command,
#                 detach=True,
#             )
#             with output_file.open("w") as outp:
#                 outp.write(" ".join(command))
#                 outp.write("\n")
#                 for line in container.logs(stream=True):
#                     line = line.strip().decode("utf-8")
#                     print(line)
#                     outp.write(line)
#                 outp.write("\n" + " ".join(command))

#         job = ppg.FileGeneratingJob(outfile, __dump).depends_on(deps)
#         job.cores_needed = 17
#         return job


def CRISPResso2(
    report_name: str,
    raw: Sample,
    df_amplicons: DataFrame,
    output_folder: Path,
    options: Optional[Dict[Any, Any]] = None,
    quantification_window_size=10,
    quantification_window_center=-3,
    dependencies: List[Job] = [],
):
    output_folder.mkdir(parents=True, exist_ok=True)
    output_file = Path(output_folder / f"CRISPResso_on_{report_name}.html")
    input_files = raw.get_aligner_input_filenames()

    def __write(output_file):
        call_crispresso2(
            report_name,
            input_files,
            df_amplicons,
            output_folder,
            options,
            quantification_window_size,
            quantification_window_center,
        )

    return (
        ppg.FileGeneratingJob(output_file, __write)
        .depends_on(raw.prepare_input())
        .depends_on(dependencies)
    )


def call_crispresso2(
    report_name: str,
    input_files: List[str],
    df_amplicons: DataFrame,
    output_folder: Path,
    options: Optional[Dict[Any, Any]] = None,
    quantification_window_size=10,
    quantification_window_center=-3,
):
    amplicon_seq = ",".join(df_amplicons.Amplicon.values)
    sgrnas = ",".join(df_amplicons.sgRNA.values)
    amplicon_names = ",".join(df_amplicons.Gen.values)
    output_folder.mkdir(parents=True, exist_ok=True)
    docker_command = [
        "docker",
        "run",
        "-v",
        "${ANYSNAKE2_PROJECT_DIR}:/project",
        "-w",
        "/project",
        "-i",
        "pinellolab/crispresso2",
    ]
    command = docker_command + [
        "CRISPResso",  # , "-h"]
        "--fastq_r1",
        str(input_files[0]),
        "--amplicon_seq",
        amplicon_seq,
        "--guide_seq",
        sgrnas,
        "--quantification_window_size",
        f"{quantification_window_size}",
        "--quantification_window_center",
        f"{quantification_window_center}",
        "--base_editor_output",
        "--exclude_bp_from_right",
        "1",
        "--exclude_bp_from_left",
        "1",
        "-o",
        str(output_folder),
        "-n",
        report_name,
        "-an",
        amplicon_names,
    ]
    # if raw_sample.is_paired:
    if len(input_files) == 2:
        command.extend(
            [
                "--fastq_r2",
                str(input_files[1]),
            ]
        )
    if options is not None:
        command.append(dict_to_string_of_items(options))
    cmd = " ".join(command)
    try:
        subprocess.run(cmd, shell=True)
    except subprocess.CalledProcessError:
        print(cmd)
        raise
