#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is a skeleton file that can serve as a starting point for a Python
console script. To run this script uncomment the following lines in the
[options.entry_points] section in setup.cfg:

    console_scripts =
         fibonacci = genome_editing.skeleton:run

Then run `python setup.py install` which will install the command `fibonacci`
inside your current environment.
Besides console scripts, the header (i.e. until _logger...) of this file can
also be used as template for Python modules.

Note: This skeleton file can be safely removed if not needed!
"""

import argparse
import sys
import logging

from genome_editing import __version__

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"

_logger = logging.getLogger(__name__)


def fib(n):
    """Fibonacci example function

    Args:
      n (int): integer

    Returns:
      int: n-th Fibonacci number
    """
    assert n > 0
    a, b = 1, 1
    for i in range(n-1):
        a, b = b, a+b
    return a

import os
import docker
import pypipegraph as ppg
from pathlib import Path


class Crispresso2:

    def __init__(self, fastq_dir = "/project/incoming"):
        self.image = "pinellolab/crispresso2"
        self.wdir = "/project"
        self.volumes = {
            os.environ["ANYSNAKE_PROJECT_PATH"]:{"bind" : "/project", "mode" :"rw"},
        }
        self.client = docker.from_env()

    def get_version(self):
        s = self.get_help()
        version = s[s.find("CRISPResso version ")+19:s.find("]")]
        return version


    def get_help(self):
        command_help = ["CRISPResso",  "-h"]
        container = self.client.containers.run(
            image = self.image, 
            volumes = self.volumes, 
            working_dir = self.wdir, 
            command = command_help
            )
        return container.decode("utf8")
    
    def print_help(self):
        print(self.get_help())

    def crispresso_run(
        self, 
        raw_sample,
        amplicon_region,
        sgRNA,
        genome,
        output_folder = None,
        quantification_window_size = 10,
        quantification_window_center = -10,
        dependencies = [],
        **options
        ):
        inputfiles = raw_sample.get_aligner_input_filenames()
        output_file = f"CRISPResso_on_{raw_sample.name}"
        deps = dependencies
        deps.append(ppg.ParameterInvariant(f"PI_{output_file}", [self.image, self.volumes, self.wdir, str(options), quantification_window_center, quantification_window_size] + amplicon_region))
        deps.append(raw_sample.prepare_input())
        outputfolder = output_folder
        if outputfolder is None:
            outputfolder = f"results/crispresso_{self.get_version()}/{raw_sample.name}"

        if isinstance(outputfolder, str):
            outputfolder = Path(outputfolder)
        outfile = outputfolder / output_file
        outfile.parent.mkdir(parents = True, exist_ok = True)
        chrm, start, stop = amplicon_region
        amp_sequence = genome.get_genome_sequence(str(chrm), start, stop)
        quantification_window_size
        def __dump():
            command = [
                "CRISPResso", 
                "--fastq_r1", 
                str(inputfiles[0]),
                "--amplicon_seq",
                amp_sequence,
                "--guide_seq",
                sgRNA, 
                "--quantification_window_size",  
                f"{quantification_window_size}", 
                "--quantification_window_center", 
                f"{quantification_window_center}", 
                "--base_editor_output",
                "-o",
                str(outputfolder),
                "-n",
                raw_sample.name
                ]
            if raw_sample.is_paired:
                    command.extend([
                    "--fastq_r2", 
                    str(inputfiles[1]),
                    ])
            for k in options:
                command.extend([k, str(options[k])])
            container = self.client.containers.run(
                self.image, 
                volumes = self.volumes, 
                working_dir = self.wdir, 
                command = command
                )
            print(container.decode("utf8"))
        return ppg.FileGeneratingJob(outfile, __dump).depends_on(deps)
