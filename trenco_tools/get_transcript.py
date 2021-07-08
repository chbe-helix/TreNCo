#!/usr/bin/env python3

'''Get TSS per transcript using Gencode annotations file for a list of biotypes.'''

import subprocess
import argparse
import os
import logging
from trenco_modules.trenco_core import process_tss, build_dir
from trenco_modules.trenco_args import txn_args

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

    txn_args(parser)

    build_dir()

    args = parser.parse_args()
    version = args.annotations
    organism = args.organism
    biotypes = args.biotypes   
    fname = args.annotfname 
    
    _, _ = process_tss(version, organism, biotypes, fname, logger)
