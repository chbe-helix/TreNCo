#!/usr/bin/env python3

'''Make matrix of log2(TPM) of enhancers for all samples.'''

import argparse
import subprocess
import os
import logging
from trenco_modules.trenco_args import merge_enh_args
from trenco_modules.trenco_core import merge_enhancer, build_dir

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

    merge_enh_args(parser)
    args = parser.parse_args()

    build_dir()

    enhancers = args.enhancers
    target = args.target
    alignment_files = args.alignments

    _, _ = merge_enhancer(enhancers, target, alignment_files, logger)
