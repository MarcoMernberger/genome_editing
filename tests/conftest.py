#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for genome_editing.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import pytest
import mbf_align
import mbf_genomes
import pypipegraph as ppg
import collections
from pypipegraph.testing.fixtures import (  # noqa:F401
    new_pipegraph,
    pytest_runtest_makereport,
)
from mbf_externals.testing.fixtures import local_store  # noqa:F401
from mbf_qualitycontrol.testing.fixtures import new_pipegraph_no_qc  # noqa:F401
from pathlib import Path
from genome_editing import PysamCheck

datapath = Path(__file__).parent / "data"

test_data = {
    "M03491:93:000000000-J8FDD:1:1107:13360:11497": {
        "R1": [
            ('1_deletion_TGAAGGAGATG_', 144, 154, -1, -1),
            ('1_mismatch__T', 238, 238, 226, 226),
        ],
        "R2": [
            ('1_deletion_TGAAGGAGATG_', 144, 154, -1, -1),
            ('1_mismatch__T', 238, 238, 226, 226),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1107:22955:11000": {
        "R1": [
            ('1_deletion_GTGATGAAGGAGATGGGA_', 140, 157, -1, -1),
            ('1_insertion__CCTC', 159, 160, 141, 144),
        ],
        "R2": [
            ('1_deletion_GTGATGAAGGAGATGGGA_', 140, 157, -1, -1),
            ('1_mismatch__G', 203, 203, -1, -1),
            ('1_mismatch__A', 205, 205, -1, -1),
            ('1_insertion__CCTC', 159, 160, 141, 144),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1105:19721:13683": {
        "R1": [
            ('1_insertion__A', 142, 143, 142, 142),
        ],
        "R2": [
            ('1_insertion__A', 142, 143, 142, 142),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1101:3803:18971": {
        "R1": [
            ('1_deletion_ATGTGATGAAGGAGATGGGAGGCCATCA_', 138, 165, -1, -1),
            ('1_insertion__G', 171, 172, 143, 143),
            ('1_insertion__AATGT', 174, 175, 147, 151),
        ],
        "R2": [
            ('1_deletion_ATGTGATGAAGGAGATGGGAGGCCATCA_', 138, 165, -1, -1),
            ('1_insertion__G', 171, 172, 143, 143),
            ('1_insertion__AATGT', 174, 175, 147, 151),
            ('1_mismatch__G', 203, 203, -1, -1),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1101:18769:19976": {
        "R1": [
            ('1_deletion_ATGTGA_', 138, 143, -1, -1),
            ('1_mismatch__A', 9, 9, -1, -1),
        ],
        "R2": [
            ('1_deletion_ATGTGA_', 138, 143, -1, -1),
        ]
    },
}


@pytest.fixture
@pytest.mark.usefixtures("new_pipegraph_no_qc")
def tgenome(tmpdir):
    ppg.new_pipegraph()
    genome = mbf_genomes.FileBasedGenome(
        "FBgenomeGTF", datapath / "genome.fasta", gtf_file=datapath / "genes.gtf", cache_dir=tmpdir
    )
    return genome


@pytest.fixture
@pytest.mark.usefixtures("new_pipegraph_no_qc")
def tlane(tgenome):
    bam_path = datapath / "test_reformat.bam"
    lane = mbf_align.lanes.AlignedSample("test_lane", bam_path, tgenome, True, "TT123")
    return lane


@pytest.fixture
def edits_to_find():
    edits = set()
    for read_name in test_data:
        count = collections.Counter()
        for edit in test_data[read_name]["R1"] + test_data[read_name]["R2"]:
            count[edit] += 1
        for edit in count:
            if count[edit] == 2:
                edits.add(edit[:3])
    return edits


@pytest.fixture
def read_edits_to_find():
    by_read = {}
    for read_name in test_data:
        count = collections.Counter()
        for edit in test_data[read_name]["R1"] + test_data[read_name]["R2"]:
            count[edit] += 1
        edits = []
        for edit in count:
            if count[edit] == 2:
                _, etype, original, replacement = edit[0].split("_")
                refstart, refstop = str(edit[1]), str(edit[2])
                qerstart, qerstop = str(edit[3]), str(edit[4])
                edits.append("_".join([etype, refstart, refstop, original, replacement, qerstart, qerstop]))
        by_read[read_name] = sorted(edits)
    return by_read


@pytest.fixture
def checknormal(tgenome):
    event_window_by_amplicon = {"1_HPRT": (30, 282), "2_TLR": (30, 317)}  # these bases are outside in 1-based igv
    primers_or_range_by_amplicon_id = {
        "2_TLR": ("agctggacggcgacgtaaac".upper(), mbf_genomes.common.reverse_complement("gatgcggttcaccagggtgt".upper())),
        "1_HPRT1": ("tctcactatattgcccaggt".upper(), mbf_genomes.common.reverse_complement("cacaatagctcttcagtctgat".upper()))
    }
    check = PysamCheck(
        primers_or_range_by_amplicon_id,
        event_window_by_amplicon,
        tgenome,
        paired=True,
        mapping_qality_filter=20,
        single_read_sufficient=False
    )
    return check
