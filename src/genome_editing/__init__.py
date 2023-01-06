# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from .crispresso import (
    CRISPResso2,
    get_amplicon_sequence,
    write_amplicon_sequences_by_name,
    call_crispresso2,
    create_crispresso_input_frame,
)
from .plots import *

# from .amplican import Amplican
# from .pycheck import PysamCheck
