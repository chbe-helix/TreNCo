import pandas as pd
import glob
import numpy as np
import os, sys
import seaborn as sb
import multiprocessing
import argparse
from trenco_modules.trenco_args import enh_gene_nw_args
from trenco_modules.trenco_core import get_enh_networks, build_dir

s_gene_tpm = None
s_enh_tpm = None
sample = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    enh_gene_nw_args(parser)

    args = parser.parse_args()

    sample = args.sample
    gene = args.gene
    enh = args.enh
    sample_out = args.output
    threads = args.threads
    enhBED = args.enhBED
    geneBED = args.geneBED
    tadBED = args.tadBED

    build_dir()

    _ = get_enh_networks(sample, 
                         gene, 
                         enh, 
                         sample_out, 
                         threads, 
                         enhBED, 
                         geneBED, 
                         tadBED)
