#!/usr/bin/env python3

'''Get universe of enhancer boundries not overlapping TSS.'''

import subprocess
import argparse
import os
import logging
from trenco_modules.trenco_args import enh_bound_args
from trenco_modules.trenco_core import enhancer_bounds, build_dir

EPILOG = '''
For more details:
        %(prog)s --help
'''

# SETTINGS

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    enh_bound_args(parser)

    args = parser.parse_args()
    tss = args.tss
    region = args.region
    sizes = args.sizes
    distance = args.distance
    peaks = args.peaks
    prange = [int(x) for x in args.promoter.split("-")]

    build_dir()

    _ = enhancer_bounds(tss, region, sizes, distance, peaks, prange, logger)
