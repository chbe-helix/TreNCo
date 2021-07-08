#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
# Wraper script for running all of trenco. The main scripts are found in      #
# trenco_core module and are imported.                                        #
#                                                                             #
# Inputs required:                                                            #
# DESIGN     tab separated file containing link information of samples        #
#               to names                                                      #
# ALIGNMENT  bam files form ChIP-seq experiments for enhancers                #
# BED        file of TAD regions to use                                       #
# EXPRESSION tsv file containing transcript expression                        #
# PEAKS bed  files with called enhancer peaks                                 #
# ANNOTAIONS gtf file vM4 recommended with gene annotations                   #
#                                                                             #
# Major Outputs:                                                              #
# Three folders:                                                              #
#   - Results - contains graph results of the networks                        #
#   - Log     - contains the logs of the scripts                              #
#   - Process - contains intermediate files and plots                         #
# EDGE csv formatted information containing graph link information between    #
#       in an array                                                           #
# NODE csv formatted information containining expression data for each gene   #
#       node                                                                  #
# --------------------------------------------------------------------------- #

import os
import json
import logging
import multiprocessing as mp
import trenco_modules.trenco_args as targs
import trenco_modules.trenco_core as tcore
from argparse import ArgumentParser

debug = False

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

_ORAGNISM_KEY = {'human' : ['Homo', 'sapiens'], 
                 'mouse' : ['Mus', 'musculus'], 
                 'worm' : ['Caenorhabditis', 'elegans'], 
                 'fly' : ['Drosophila', 'melanogaster']}

_ORAGNISM_GENOME = {'human' : 'hg38', 
                    'mouse' : 'mm10'}

_GENOME_CONVERSION = {'hg38' : 'GRCh38.p13',
                      'mm10' : 'GRCm38.p6'}

def initial_setup(version,
                  organism,
                  biotypes,
                  annotfile,
                  region, 
                  distance, 
                  peaks, 
                  prange,
                  expression_files,
                  alignment_files,
                  target,
                  gvers,
                  meme_db,
                  memedb,
                  logger):

    ### Build Directories for analysis if they don't exsist
    print("Building Core Directories")
    tcore.build_dir()

    ### Get Transcripts
    print("Gathering Transcripts")
    tss, bedgene = tcore.process_tss(version,
                                     organism,
                                     biotypes,
                                     annotfile,
                                     logger)
    
    ### Bounding Enhancers
    print("Bounding Enhancers")
    mergef, probed = tcore.enhancer_bounds(tss, 
                                           region, 
                                           gvers, 
                                           distance, 
                                           peaks, 
                                           prange, 
                                           logger)

    ### Building transcript log2 matrix
    print("Building Transcript Matrix")
    txn_log2TPM = tcore.merge_transcript(expression_files,
                                         logger)

    ### Building enhancer log2 matrix
    print("Building Enhancer Matrix")
    enh_log2TPM, bedenh = tcore.merge_enhancer(mergef, 
                                               target, 
                                               alignment_files, 
                                               logger)

    if organism in _ORAGNISM_KEY:
        organism = _ORAGNISM_KEY[organism]
        if not memedb:
            memedb = organism

    ### Building transcript log2 matrix
    print("Building Transcription Factor Matricies")
    mtcs = []
    beds = [probed, bedenh]

    if gvers in _GENOME_CONVERSION:
        gvers = _GENOME_CONVERSION[gvers]

    for bed in beds:
        mxfn = tcore.build_tf_matrix(bed, 
                                    meme_db, 
                                    memedb, 
                                    " ".join(organism),
                                    gvers, 
                                    logger)
        mtcs.append(mxfn)

    return (bedgene, bedenh, txn_log2TPM, enh_log2TPM, mtcs[0], mtcs[1])

def _DEBUG():
    bedgene     = "process/gene.bed" 
    bedenh      = "process/enh.bed" 
    txn_log2TPM = "process/transcript_expression_log2TPM_matrix.txt" 
    enh_log2TPM = "process/enhancer_H3K27ac_log2TPM_signal_matrix.txt" 
    mtx1        = "process/TFxgencode.vM4.annotation_capped_sites_"\
                        "1000u200d_promoter.bed_score_matrix.txt" 
    mtx2        = "process/TFxenh.bed_score_matrix.txt"
    return (bedgene, bedenh, txn_log2TPM, enh_log2TPM, mtx1, mtx2)

def full_build(samples,
               gene,
               enh,
               enhBED,
               geneBED,
               tadBED,
               enh_mtx,
               gene_mtx,
               threads):
    for sample, mask in samples.items():
        gene_vect = tcore.get_enh_networks(sample, 
                                           gene, 
                                           enh, 
                                           sample, 
                                           mp.cpu_count()-1, 
                                           enhBED, 
                                           geneBED, 
                                           tadBED,
                                           mask)

        tcore.tis_enh_gene_net(sample,
                               sample,
                               enh_mtx,
                               gene_mtx,
                               gene_vect,
                               threads)

if __name__ == '__main__':
    parser = ArgumentParser(
        description = "TAD aware Regulatory Network Construction (TReNCo)",
        usage="trenco --design [DESIGN FILE (txt)] \
                --alignment [ALIGNMENT FILES (tsv)] \
                --expression [EXPRESSION FILES (bam)] \
                --peaks [PEAK FILES (bed)] -g [GENOME VERSION] [OPTIONS]")

    targs.full_trenco_args(parser)
    targs.enh_bound_args(parser, False)
    targs.txn_args(parser)
    targs.tf_mtx_args(parser, False)
    
    args = parser.parse_args()

    files = [args.peaks, 
             args.alignment, 
             args.expression, 
             args.tadBED, 
             args.design]
    for fi in files:
        if not isinstance(fi, list):
            fi = [fi]

        for f in fi:
            if not os.path.exists(f):
                logger.warning("Cannot find {}".format(f))
                exit(1)
    
    peaks, alignments, expressions, tads, design = files
    prange = [int(x) for x in args.promoter_range.split("-")]

    design_key     = {}
    design_key_inv = {}
    with open(design, 'r') as ifi:
        for line in ifi:
            line = line.split()
            design_key[line[0]] = line[1:]
            for l in line[1:]:
                design_key_inv[l] = line[0]

    if args.db:
        memedb = args.db
    else:
        memedb = ''

    if not args.gvers:
        if args.organism not in _ORAGNISM_GENOME:
            print("Please specify a genome version for this species")
            exit(0)
        args.gvers = _ORAGNISM_GENOME[args.organism]

    if not debug:
        geneBED, \
          enhBED, \
          txn_log2TPM, \
          enh_log2TPM, \
          tf_gene_mx, \
          tf_enh_mx \
            = initial_setup(args.annotations,
                            args.organism,
                            args.biotypes,
                            args.annotfname,
                            args.region, 
                            args.distance, 
                            peaks, 
                            prange,
                            args.expression,
                            args.alignment,
                            args.target,
                            args.gvers,
                            args.meme_db,
                            memedb,
                            logger)
    else:
        geneBED, \
          enhBED, \
          txn_log2TPM, \
          enh_log2TPM, \
          tf_gene_mx, \
          tf_enh_mx \
            = _DEBUG()    
    
    full_build(design_key,
               txn_log2TPM,
               enh_log2TPM,
               enhBED,
               geneBED,
               args.tadBED,
               tf_enh_mx,
               tf_gene_mx,
               mp.cpu_count()-1)