#!/usr/bin/env python3

'''Make matrix of log2(TPM) of transcripts for all samples.'''

import argparse
import os
import logging
from trenco_modules.trenco_args import merge_txn_args
from trenco_modules.trenco_core import merge_transcript, build_dir

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

    merge_txn_args(parser)
    args = parser.parse_args()

    build_dir()

    expression_files = args.expression
    _ = merge_transcript(expression_files, logger)
