#!/usr/bin/env python3

import argparse, os

"""
Trenco Module for arguments
"""

def txn_args(parser):
    parser.add_argument('--annotation-file',
                        dest = 'annotfname',
                        default = '',
                        help="Genode annotations file in gtf format (overwrites --annotation-version and --organism")

    parser.add_argument('--annotation-version',
                        dest="annotations",
                        default="vM4",
                        help="The Gencode annotations file in gtf format. (Default: vM4) (WARNING: All entries are indexed to this version)")

    parser.add_argument('--organism',
                        default="mouse",
                        help="Organism gencode to download (Default: mouse)"
                        )

    parser.add_argument('-b', '--biotypes',
                        help="The biotypes to get transcript TSS. (default: protein)",
                        nargs='+',
                        default=['protein_coding'])

def enh_bound_args(parser, tot = True):
    if tot:
        parser.add_argument('-t', '--tss',
                            help="The Gencode TSS file.",
                            required=True)
        parser.add_argument('-s', '--sizes',
                            help="The chrome sizes file or genome number (ie mm10)",
                            required=True)
    parser.add_argument('-p', '--peaks',
                        help="The full path to peak files in bed format",
                        nargs='+',
                        required=True)

    #parser.add_argument('--geneGTF',
    #                    help="GTF file of genes from Gencode (Default gene_txn.gtf from get_trancript script)",
    #                    default = "gene_txn.gtf") 

    parser.add_argument('-r', '--region',
                        help="The number of bases pairs to exclude around TSS (Default: 2500)",
                        type=int,
                        default=2500)

    parser.add_argument('-q', '--promoter-range',
                        help="Range of numbers before TSS and after TSS to consider as Promoter (Default: 1000-200)",
                        type=str,
                        default="1000-200")

    parser.add_argument('-d', '--distance',
                        help="The number of bases pairs to merge between adjacent peaks (Default: 150)",
                        type=int,
                        default=150)

def merge_txn_args(parser):
    parser.add_argument('-e', '--expression',
                        help="The full path to peak files in tsv format",
                        nargs='+',
                        required=True)

def merge_enh_args(parser):
    parser.add_argument('-e', '--enhancers',
                        help="The universe of enhancer files.",
                        required=True)

    parser.add_argument("-t", "--enhMarks",
                        dest='target',
                        type=str,
                        default="H3K27ac",
                        help="Mark for enchancers: Default H3K27ac")

    parser.add_argument('-a', '--alignments',
                        help="The full path to sorted alignment files in bam format.",
                        nargs='+',
                        required=True)

def full_trenco_args(parser):
    path = '/'.join(os.path.realpath(__file__).split('/')[:-2])
    parser.add_argument("--design",
                        type=str,
                        required=True,
                        help="Design file containing link information to samples.")
    parser.add_argument("--alignment",
                        nargs='+',
                        required=True,
                        help="Full path to ChIP alingment files in bam format")
    parser.add_argument("--expression",
                        nargs='+',
                        required=True,
                        help="Full path to transcript expression table in tsv format")
    parser.add_argument("--enhMarks",
                        dest='target',
                        type=str,
                        default="H3K27ac",
                        help="Mark for enchancers: Default H3K27ac")    
    parser.add_argument("--tadBED",
                        type=str,
                        default="%s/common_files/TAD_positions_mm10_sorted.bed" % path,
                        help="TAD file: Default - mm10 TAD in common_files")


def tf_mtx_args(parser, spec = True):
    parser.add_argument("--meme-db",
                        type=str,
                        default="cis-bp",
                        help="MEME database to use (Default: cis-bp)")
    parser.add_argument("--db",
                        type=str,
                        help="Motif database name if different from species (ex JASPER CORE 2014 for jasper)")
    if spec:
        parser.add_argument('-s', '--species',
                            dest='refID',
                            nargs='+',
                            required=True,
                            help = "Scientific name of organism (can use these names for ease: human, mouse, fly, worm)")
    parser.add_argument('-g', '--genome-version',
                        dest='gvers',
                        type=str,
                        help = "Version of genome to use. Default is newest")
    parser.add_argument('--bed',
                        dest='bed',
                        type=str,
                        help = "ChIP and Promoter bed file for getting motifs (ex enh.bed,promoter.bed)")

def enh_gene_nw_args(parser):
    parser.add_argument("-e", "--enh", help="enhancer by samples log2 TPM quantification matrix", type=str)
    parser.add_argument("-g", "--gene", help="gene by samples log2 TPM quantification matrix", type=str)
    parser.add_argument("-ta", "--tadBED", help='sorted tad annotation in bed file format', type=str)
    parser.add_argument("-ga", "--geneBED", help='gene annotation in bed file format', type=str)
    parser.add_argument("-ea", "--enhBED", help='enh annotation in bed file format')
    parser.add_argument("-s", "--sample", help='sample to construct the network', type=str)
    parser.add_argument("-o", "--output", help="output directory", type=str)
    parser.add_argument("-p", "--threads", help="Threads", type=int, default=30)

def tis_gene_networks(parser):
    parser.add_argument("-d", "--dir", help="directory containing the output of get_enh_gene_networks.py", type=str)
    parser.add_argument("-s", "--sample", help='sample to construct the network', type=str)
    parser.add_argument("-p", "--threads", help='Number of threads to use', type=int, default=30)
    parser.add_argument("-x1", "--matrix1", help='TF by enchancer matrix file path', type=str)
    parser.add_argument("-x2", "--matrix2", help="TF by gene promoter matrix file path", type=str)
    parser.add_argument("-v", "--vector", help="Expression vector for the sample from RNA-seq", type=str)