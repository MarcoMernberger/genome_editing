#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from genome_editing.crispresso import Crispresso2

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def test_print():
    crip = Crispresso2()
    help_text = crip.get_help()
    assert "usage: CRISPResso" in help_text



