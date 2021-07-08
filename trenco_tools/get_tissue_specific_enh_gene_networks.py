import multiprocessing
import argparse
from trenco_modules.trenco_args import tis_gene_networks
from trenco_modules.trenco_core import tis_enh_gene_net, build_dir

sample = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    tis_gene_networks(parser)

    args = parser.parse_args()

    sample = args.sample
    threads = args.threads
    outdir = args.outdir
    mx1 = args.matrix1
    mx2 = args.matrix2
    vect = args.vector

    build_dir()

    tis_enh_gene_net(sample,
                     outdir,
                     mx1,
                     mx2,
                     vect,
                     threads)
