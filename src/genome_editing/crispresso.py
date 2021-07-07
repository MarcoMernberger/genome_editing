#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import random
import docker
import pypipegraph as ppg
import mbf_align
import mbf_genomes
from pathlib import Path
from typing import List, Any


class Crispresso2:
    """
     Wrapper for the crispresso2 docker.
    
    [extended_summary]
    
    Returns
    -------
    [type]
        [description]
    """

    def __init__(self, fastq_dir="/project/incoming"):
        self.image = "pinellolab/crispresso2"
        self.wdir = "/project"
        self.volumes = {
            os.environ["ANYSNAKE_PROJECT_PATH"]: {"bind": "/project", "mode": "rw"},
        }
        self.client = docker.from_env()

    def get_version(self):
        s = self.get_help()
        version = s[s.find("CRISPResso version ") + 19 : s.find("]")]
        return version

    def get_help(self):
        command_help = ["CRISPResso", "-h"]
        container = self.client.containers.run(
            image=self.image,
            volumes=self.volumes,
            working_dir=self.wdir,
            command=command_help,
            detach=True,
        )
        ret = ""
        for line in container.logs(stream=True):
            ret += line.decode("utf-8")

        return ret

    def print_help(self):
        print(self.get_help())

    def crispresso_run(
        self,
        raw_sample: mbf_align.raw.Sample,
        amplicon_regions: List[List[Any]],
        sgrnas: List[str],
        genome,
        names: List[str] = None,
        output_folder=None,
        quantification_window_size=10,
        quantification_window_center=-3,
        dependencies=[],
        options={},
    ):
        inputfiles = raw_sample.get_aligner_input_filenames()
        if names is None:
            names = []
            for sg in sgrnas:
                names.append(f"{raw_sample.name.replace('.', '_')}_{sg}")
        outputfolder = output_folder
        if outputfolder is None:
            outputfolder = f"results/crispresso_{self.get_version()}/{raw_sample.name}"
        if isinstance(outputfolder, str):
            outputfolder = Path(outputfolder)
        aname = ",".join(names)
        name = aname.replace(",", "_")
        output_file = f"CRISPResso_on_{name}"
        outfile = outputfolder / output_file
        outfile.parent.mkdir(parents=True, exist_ok=True)
        amp_sequences = []
        for amplicon_region in amplicon_regions:
            chrm, start, stop = amplicon_region
            amp_sequence = genome.get_genome_sequence(str(chrm), start, stop)
            amp_sequences.append(amp_sequence)
        amp_sequence = ",".join(amp_sequences)
        sgrna = ",".join(sgrnas)
        deps = dependencies
        deps.append(
            ppg.ParameterInvariant(
                f"PI_{outfile}",
                [
                    self.image,
                    self.volumes,
                    self.wdir,
                    str(options),
                    quantification_window_center,
                    quantification_window_size,
                ]
                + amplicon_regions
                + sgrnas
                + names,
            )
        )
        deps.append(raw_sample.prepare_input())

        def __dump():
            command = [
                "CRISPResso",
                "--fastq_r1",
                str(inputfiles[0]),
                "--amplicon_seq",
                amp_sequence,
                "--guide_seq",
                sgrna,
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
                str(outputfolder),
                "-n",
                name,
                "-an",
                aname,
            ]
            if raw_sample.is_paired:
                command.extend(
                    ["--fastq_r2", str(inputfiles[1]),]
                )
            for k in options:
                if len(str(options[k])) == 0:
                    command.append(k)
                else:
                    command.extend([k, str(options[k])])
            time.sleep(random.random() * 2)
            container = self.client.containers.run(
                self.image,
                volumes=self.volumes,
                working_dir=self.wdir,
                command=command,
                detach=True,
            )
            with open(outputfolder / (output_file + ".pylog"), "w") as outp:
                for line in container.logs(stream=True):
                    line = line.strip().decode("utf-8")
                    print(line)
                    outp.write(line)
                outp.write("\n" + " ".join(command))

        job = ppg.FileGeneratingJob(outfile, __dump).depends_on(deps)
        job.cores_needed = 17
        return job
