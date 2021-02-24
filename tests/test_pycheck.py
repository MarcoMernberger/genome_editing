#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import mbf_align
import pypipegraph as ppg
import mbf_genomes
import pandas as pd
import pysam
from pathlib import Path
from genome_editing import PysamCheck

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_lane(tlane, tgenome):
    assert tlane.name == "test_lane"


def test_init(checknormal, tgenome):
    assert checknormal.name == "PysamCheck"
    assert checknormal.genome.name == tgenome.name
    assert checknormal.paired
    assert checknormal.discard_if_pairs_differ
    assert checknormal.types_to_consider == ["insertion", "deletion", "mismatch"]


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_normal(tlane, tgenome, tmpdir, edits_to_find, checknormal, read_edits_to_find):
    checknormal.count_events(tlane, Path(tmpdir))
    ppg.run_pipegraph()
    df_per_event = pd.read_csv(tmpdir / f"Pycheck_{tlane.name}_per_event.tsv", sep="\t")
    df_per_event.fillna("", inplace=True)
    df_per_read = pd.read_csv(tmpdir / f"Pycheck_{tlane.name}_per_read.tsv", sep="\t")
    df_tallied = pd.read_csv(tmpdir / f"Pycheck_{tlane.name}_tallied_reads.tsv", sep="\t")
    # test per event results
    for _, row in df_per_event.iterrows():
        edit = ("_".join([str(row["Counts"]), row["Type"], str(row["Original"]), str(row["Replacement"])]), row["Start"], row["Stop"])
        assert edit in edits_to_find
    # test per read results
    assert len(df_per_read) == len(read_edits_to_find)
    for _, row in df_per_read.iterrows():
        edits_found = sorted(row["Info"].split(","))
        if len(edits_found) > 1:
            etype = "multi"
        else:
            etype = edits_found[0].split("_")[0]
        assert row["Type"] == etype
        assert edits_found in read_edits_to_find.values()
    # test tallied
    assert len(df_tallied) == len(read_edits_to_find)
    for _, row in df_tallied.iterrows():
        assert len(row["Sequence"].split(",")) == 2
        assert row["Count"] == 1


def test_dependencies(checknormal):
    deps = checknormal.get_dependencies()
    for x in deps:
        print(type(x))
        assert (isinstance(x, ppg.job.Job) or isinstance(x, list))


def test_tally_reads(checknormal, tlane):
    counter, read_dict = checknormal.tally_reads(tlane)
    assert len(counter) == 5
    for sequences in counter:
        assert isinstance(sequences[0], str)
        assert isinstance(sequences[1], str)
        assert counter[sequences] == 1
    for aligned_segment_tuple in read_dict.values():
        assert isinstance(aligned_segment_tuple[0], pysam.AlignedSegment)
        assert isinstance(aligned_segment_tuple[1], pysam.AlignedSegment)


def test_get_homology(checknormal):
    deletion = "ATCGCT"
    flanking_forward = "ATGGGG"
    flanking_reverse = "GGGCT"
    hom1 = checknormal._get_homology(deletion, flanking_forward, "")
    hom2 = checknormal._get_homology(deletion, "", flanking_reverse)
    hom3 = checknormal._get_homology(deletion, flanking_forward, flanking_reverse)
    assert hom1 == "AT"
    assert hom2 == "GCT"
    assert hom3 == hom2

"""
def test_eval(checknormal):
    cigars = ['9S128=28D6=1I3=5I99=', '95=28D6=1I3=5I28=1X100=11S'] 
    aln_pos_start = [11, 44] 
    query_seq = [
        'TCACTGGCATCTCACTATATTGCCCAGGTTGGTGTGGAAGTTTAATGACTAAGAGGTGTTTGTTATAAAGTTTAATGTATGAAACTTTCTATTAAATTCCTGATTTTATTTCTGTAGGACTGAACGTCTTGCTCGAGCATTGTGAGCAATGTCCTCTGTGTGCTCAAGGGGGGCTATAAATTCTTTGCTGACCTGCTGGATTACATCAAAGCACTGAATAGAAATAGTGATAGATCCATTCCTATGACTGT', 
        'TAATGACTAAGAGGTGTTTGTTATAAAGTTTAATGTATGAAACTTTCTATTAAATTCCTGATTTTATTTCTGTAGGACTGAACGTCTTGCTCGAGCATTGTGAGCAATGTCCTCTGTGTGCTCAAGGGGGGCTATAAAGTCTTTGCTGACCTGCTGGATTACATCAAAGCACTGAATAGAAATAGTGATAGATCCATTCCTATGACTGTAGATTTTATCAGACTGAAGAGCTATTGTGTGCCAGCAGCTA'
        ]
    ref = [
        'NNNNNNNNNNTCTCACTATATTGCCCAGGTTGGTGTGGAAGTTTAATGACTAAGAGGTGTTTGTTATAAAGTTTAATGTATGAAACTTTCTATTAAATTCCTGATTTTATTTCTGTAGGACTGAACGTCTTGCTCGAGATGTGATGAAGGAGATGGGAGGCCATCACATTGTAGCCCTCTGTGTGCTCAAGGGGGGCTATAAATTCTTTGCTGACCTGCTGGATTACATCAAAGCACTGAATAGAAATAGTGATAGATCCATTCCTATGACTGTAGATTTTATCAGACTGAAGAGCTATTGTGNNNNNNNNNN',
        'NNNNNNNNNNTCTCACTATATTGCCCAGGTTGGTGTGGAAGTTTAATGACTAAGAGGTGTTTGTTATAAAGTTTAATGTATGAAACTTTCTATTAAATTCCTGATTTTATTTCTGTAGGACTGAACGTCTTGCTCGAGATGTGATGAAGGAGATGGGAGGCCATCACATTGTAGCCCTCTGTGTGCTCAAGGGGGGCTATAAATTCTTTGCTGACCTGCTGGATTACATCAAAGCACTGAATAGAAATAGTGATAGATCCATTCCTATGACTGTAGATTTTATCAGACTGAAGAGCTATTGTGNNNNNNNNNN'
        ]
    read_id = ['R1', 'R2']
    mapq = [60, 60] 
    seqnames = ['1_HPRT1', '1_HPRT1']
    strands = ['+', '-']
    counts = [1, 1]
    events = checknormal.evaluate_cigar(
        cigars,
        aln_pos_start,
        query_seq,
        ref,
        read_id,
        mapq,
        seqnames,
        strands,
        counts,
    )
    expected = [
        1_HPRT1    166  166  ...   mismatch      R2    60      1
        1_HPRT1    138  165  ...   deletion      R1    60      1
        1_HPRT1    138  165  ...   deletion      R2    60      1
        1_HPRT1    162  162  ...  insertion      R1    60      1
        1_HPRT1    166  170  ...  insertion      R1    60      1
        1_HPRT1    129  129  ...  insertion      R2    60      1
        1_HPRT1    133  137  ...  insertion      R2    60      1
    ]
    deletions = cigars[0].count("D") + cigars[1].count("D")
    insertions = cigars[0].count("I") + cigars[1].count("I")
    print(events)
    for _, row in events.iterrows():s
        assert index seqnames  start  end  ...       type read_id score counts
    raise ValueError()
"""