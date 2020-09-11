#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


test_data = {
    "M03491:93:000000000-J8FDD:1:1107:13360:11497": {
        "R1": [
            ('1_deletion_TGAAGGAGATG_', 144, 155),
            ('1_mismatch__T', 238, 238),
        ],
        "R2": [
            ('1_deletion_TGAAGGAGATG_', 144, 155),
            ('1_mismatch__T', 238, 238),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1107:22955:11000": {
        "R1": [
            ('1_deletion_GTGATGAAGGAGATGGGA_', 140, 157),
            ('1_insertion__CCTC', 159, 160),
        ],
        "R2": [
            ('1_deletion_GTGATGAAGGAGATGGGA_', 140, 157),
            ('1_mismatch__G', 203, 203),
            ('1_mismatch__A', 205, 205),
            ('1_insertion__CCTC', 159, 160),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1105:19721:13683": {
        "R1": [
            ('1_insertion__A', 142, 143),
        ],
        "R2": [
            ('1_insertion__A', 142, 143),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1101:3803:18971": {
        "R1": [
            ('1_deletion_ATGTGATGAAGGAGATGGGAGGCCATCA_', 138, 165),
            ('1_insertion__G', 171, 172),
            ('1_insertion__AATGT', 174, 175),
        ],
        "R2": [
            ('1_deletion_ATGTGATGAAGGAGATGGGAGGCCATCA_', 138, 165),
            ('1_insertion__G', 171, 172),
            ('1_insertion__AATGT', 174, 175),
            ('1_mismatch__G', 203, 203),
        ]
    },
    "M03491:93:000000000-J8FDD:1:1101:18769:19976": {
        "R1": [
            ('1_deletion_ATGTGA_', 138, 143),
            ('1_mismatch__A', 9, 9),
        ],
        "R2": [
            ('1_deletion_ATGTGA_', 138, 143),
        ]
    },
}

