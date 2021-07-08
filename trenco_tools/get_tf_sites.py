#!/usr/bin/env python3

import os, subprocess
from argparse import ArgumentParser
import logging
from trenco_modules.trenco_core import build_tf_matrix, build_dir
from trenco_modules.trenco_args import tf_mtx_args

EPILOG = '''
For more details:
        %(prog)s --help
'''

_ORAGNISM_KEY = {'human' : ['Homo', 'sapiens'], 
                 'mouse' : ['Mus', 'musculus'], 
                 'worm' : ['Caenorhabditis', 'elegans'], 
                 'fly' : ['Drosophila', 'melanogaster']}

# SETTINGS
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

if __name__ == '__main__':
    parser = ArgumentParser(
    description = "Generate TF by Gene and TF by Enhancer matricies")

    tf_mtx_args(parser)
    args = parser.parse_args()

    if args.refID[0] in _ORAGNISM_KEY:
        orgname = _ORAGNISM_KEY[args.refID[0]]
    else:
        orgname = args.refID

    if args.db:
        memedb = args.db
    else:
        memedb = orgname

    if not args.bed or ',' not in args.bed:
        print("Need two bed files, promoter and enh")
        exit(1)

    beds = args.bed.split(',')
    meme_db = args.meme_db
    gvers = args.gvers

    _ = build_tf_matrix(beds, meme_db, memedb, orgname, gvers, logger)
