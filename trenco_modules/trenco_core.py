#!/usr/bin/env python3

"""
Trenco Modules for running scripts
"""

import sys
import re
import os
import glob
import subprocess
import json
import logging
import time
import shlex
import pandas as pd
import numpy as np
import seaborn as sb
import multiprocessing as mp
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# Switch
debug     = True
tf_enh_mx = None

# --------------------------------------------------------------------------- #
# General scripts that are used in trenco wrapper and in multiple places      #
#   throughout the core module                                                #
# --------------------------------------------------------------------------- #
def build_dir():
    dirs = ['log', 'process', 'results']
    for d in dirs:
        if not os.path.exists(d):
            print("Building {}".format(d))
            os.mkdir(d)

def plot_mtx(fnmtx):
    mtx = pd.read_csv(fnmtx, sep='\t')

    try:
        g = sb.pairplot(mtx,diag_kind='kde')
    except:
        mtx.columns = mtx.columns.str.split("_").str[-2]
        g = sb.pairplot(mtx, diag_kind='kde')
        pass

    for ax in g.axes.flat:
        ax.set_xlabel(xlabel = ax.xaxis.get_label_text(), rotation=90)
        ax.set_ylabel(ylabel = ax.yaxis.get_label_text(), rotation=80)

    g.savefig("results/{}.png".format(fnmtx.split('/')[-1]))
    plt.clf()

def plot_corr_heatmap(df, df_name, size):
    """ Plot heatmap of numeric pandas DataFrame """
    a4_dims = (size, size)
    _, ax = plt.subplots(figsize=a4_dims)
    
    sb.set_context('paper', font_scale=1.8)
    g = sb.heatmap(df, yticklabels=False, cmap='Blues', annot=True, ax=ax)
    fig = g.get_figure()
    fig.savefig('{}-heatmap.png'.format(df_name))
    plt.clf()

def files_exist(fnames = []):
    for name in fnames:
        if not os.path.exists(name):
            return False

    return True

def check_if_commons(commonfile, exact = True):
    filepath = "/".join(os.path.realpath(__file__).split("/")[:-2]) \
                    + "/common_files/*"
    commons = glob.glob(filepath)
   
    for file_full in commons:
        file_name = file_full.split("/")[-1]
        if exact:
            if commonfile == file_name:
                return file_full
        else:
            if commonfile in file_name:
                return file_full

    return None

def load_convert_name(fn):
    converter = {}
    with open(fn, "r") as ifi:
        for line in ifi:
            line = line.strip()
            if line.startswith("#"):
                continue

            key, val = line.split()
            converter[key] = val

    return converter

def load_gtf_conversion(gtf, conv_from, conv_to):
    if not os.path.exists(gtf):
        print("File %s cannot be found" % gtf)
        exit(1)

    converter = {}
    with open(gtf, "r") as ifi:
        for line in ifi:
            line = line.strip()
            if line.startswith("#"):
                continue
            
            keys = line.split("\t")[-1]
            keys = [k.strip().replace("\"", "") for k in keys.split(";")]

            keyset = {}
            for k in keys:
                if not k:
                    continue

                flag, entry = k.split()
                keyset[flag] = entry

            if keyset[conv_from] not in converter:
                converter[keyset[conv_from]] = keyset[conv_to]
            else:
                if converter[keyset[conv_from]] == keyset[conv_to]:
                    continue
                elif isinstance(converter[keyset[conv_from]], list):
                    converter[keyset[conv_from]].append(keyset[conv_to])
                else:
                    converter[keyset[conv_from]] = [converter[keyset[conv_from]], 
                                                    keyset[conv_to]]

    return converter

# --------------------------------------------------------------------------- #
# Methods for Getting Transcripts from ANNOTATIONS gtf and requires           #
#   AWK script                                                                #
# --------------------------------------------------------------------------- #
def get_tss(annotation, biotypes, logger):
    '''Get the TSS.'''
    # Make biotypes file
    biotypes_file = open('process/biotypes.txt', 'w')
    biotypes_file.write("\n".join(biotypes))
    biotypes_file.close()

    # Get file basename
    file_name = os.path.basename(annotation)
    basename  = file_name.rpartition('.gtf')[0]

    # Get the gencode
    rel_path         = os.path.dirname(os.path.realpath(__file__))
    make_gencode_tss = '{}/../trenco_tools'.format(rel_path)
    tss_program      = 'make_TSS_file_from_annotation_with_confidence_better.sh'
    tss_command      = '{}/{} {} process/biotypes.txt'.format(make_gencode_tss, 
                                                              tss_program, 
                                                              annotation)

    logger.info("Running TSS caluclation %s", tss_command)

    tss = subprocess.Popen(tss_command, shell=True)
    out, err = tss.communicate()

    # Remove extra files
    os.remove('{}_exons_most5p.gff'.format(basename))
    os.remove('{}_capped_sites_nr_with_confidence.gff'.format(basename))

def process_tss(version, organism, biotypes, fname, logger):
    ''' Process the TSS into a bed file that can be used '''
    # If files don't exists
    if fname:
        annotations = fname
        if not os.path.exists(annotations):
            print("ERROR: %s not found" % annotations,
                  file=sys.stderr)
            exit(1)

    else:
        repo_fn     = "gencode.%s.annotation.gtf" % (version)
        annotations = "process/" + repo_fn
        if not os.path.exists(annotations):
            cmds = ["wget ftp://ftp.ebi.ac.uk/pub/databases/"\
                        "gencode/Gencode_%s/release_%s/%s.gz -P process" 
                            % (organism, version.replace('v',''), repo_fn),
                    "gunzip %s.gz" % (annotations)]
            
            for cmd in cmds:
                #os.system(cmd) # the code needs to wait at this point
                subprocess.check_call(cmd, shell=True)
            # this is a temporary fix to make sure the gunzip works and has
            # time to finish
            time.sleep(30)
            
    if annotations.endswith('.gtf'):
        tail = '.gtf'
    else:
        tail = '.gff'
                
    cs_annot = annotations.replace(tail, '_capped_sites.gff').split('/')[-1]    
    ofn = "process/gene.bed"

    if files_exist(['process/'+cs_annot, ofn]):
        return 'process/'+cs_annot, ofn

    # Generate Converstion table for gene names
    convertion = load_gtf_conversion(annotations, "gene_id", "gene_name")
    
    with open("process/converter.json", "w") as ofo:
        json.dump(convertion, ofo, indent = 2)

    # Create a file handler
    handler = logging.FileHandler('log/get_transcript.log')
    logger.addHandler(handler)

    # Extract biotype genes
    if os.path.exists("process/gene_txn.gtf"):
        os.remove("process/gene_txn.gtf")

    for biotype in biotypes:
        grep_cmd = "grep -w gene {} | grep {} "\
                        ">> process/gene_txn.gtf".format(annotations, biotype)

        print(grep_cmd)
        subprocess.check_call(grep_cmd, shell=True)

    subprocess.check_call("gtf2bed < process/gene_txn.gtf "\
                                "> \process/gene_txn.bed", shell=True)

    ofo = open(ofn, "w")
    with open("process/gene_txn.bed", "r") as ifn:
        for line in ifn:
            line     = line.strip().split('\t')
            idx      = line[-1]
            idx      = idx.split(';')[1].split("\"")[1]
            nline    = line[0:6]
            nline[3] = idx
            ofo.write("\t".join(nline) + "\n")

    ofo.close()
    os.system("rm process/gene_txn.gtf")
    os.system("rm process/gene_txn.bed")

    # Run get_tss
    get_tss(annotations, biotypes, logger)
    assert os.path.exists(cs_annot), 'cannot find %s' % cs_annot

    os.system('mv %s process/' % cs_annot)
    return 'process/'+cs_annot, ofn

# --------------------------------------------------------------------------- #
# Method for making Enhancer Boundries from PEAKS files and ensuring no       #
#   overlaps with TSS                                                         #
# --------------------------------------------------------------------------- #
def extend_tss(annotation, 
               sizes, 
               region, 
               prange, 
               logger):
    '''Extend TSS for masking.'''
    file_name = os.path.basename(annotation)
    if 'gtf' in annotation:
        basename = file_name.rpartition('.gtf')[0]
    else:
        if 'gff' not in annotation:
            print("Please provide a GTF or GFF file for -t option",
                  file=sys.stderr)
            exit(1)
        basename = file_name.rpartition('.gff')[0]
 
    # outputs
    up, down = prange
    promoter_file = 'process/{}_{}u{}d_promoter.bed'.format(basename, up, down)
    flanking_file = 'process/{}_{}.bed'.format(basename, region)

    if files_exist([promoter_file, flanking_file]):
        return flanking_file, promoter_file

    # Convert gff to bed file
    bed_file = 'process/{}.bed'.format(basename)
    extract = False
    if not os.path.exists(bed_file):
        extract = True
        if 'gtf' in annotation:
            convert_command = 'gtf2bed < %s > %s' % (annotation, bed_file)
        else:
            convert_command = 'gff2bed < %s | awk \'BEGIN{FS=OFS="\t"} '\
                                    '{print $1,$2,$3,$5,$4,$6}\' > %s' \
                                        % (annotation, bed_file)
    else:
        convert_command = "File Exists"    

    logger.info("Running convert with %s", convert_command)
    if extract:
        #subprocess.check_call(shlex.split(convert_command))
        os.system(convert_command)

    # Extend bed file
    flanking_command = 'bedtools slop -i {} -g {} -b {}'.format(bed_file, 
                                                                sizes, 
                                                                region)

    logger.info("Running flanking around TSS with %s", flanking_command)
    subprocess.check_call(shlex.split(flanking_command), 
                          stdout=open(flanking_file, 'w'))
    
    #Generate promoter files
    promoter_command = 'bedtools slop -i {} -g {} -l {} -r {} -s'.format(bed_file,
                                                                         sizes, 
                                                                         up, 
                                                                         down)

    logger.info("Running promoter around TSS with %s", promoter_command)
    subprocess.check_call(shlex.split(promoter_command), 
                          stdout=open(promoter_file, 'w'))

    os.remove(bed_file)
    return flanking_file, promoter_file


def merge_peaks(peaks, logger, universe=None):
    '''Merge peak files.'''
    tmp_file = 'process/tmp.bed'
    new_universe = 'process/universe_peaks.bed'

    # Merge peaks files and master file if present
    if 'gz' not in peaks:
        if universe is not None:
            merge_command = 'bedops --everything {} {} '\
                                '> {}'.format(peaks, 
                                              universe, 
                                              tmp_file)
        else:
            merge_command = 'bedops --everything {} '\
                                '> {}'.format(peaks, 
                                            tmp_file)
    else:
        if universe is not None:
            merge_command = 'gunzip -c {} '\
                                '| bedops --everything - {} '\
                                    '> {}'.format(peaks,universe,
                                                  tmp_file)
        else:
            merge_command = 'gunzip -c {} '\
                                '| bedops --everything - '\
                                    '> {}'.format(peaks, 
                                                  tmp_file)
 
    logger.info("Running merge with %s", merge_command)
    subprocess.check_call(merge_command, shell=True)
    os.rename(tmp_file, new_universe)
    return new_universe

def exclude_tss(peaks, 
                tss, 
                distance, 
                logger):
    '''Exclude peaks that overlap tss.'''
    tmp_file    = 'process/tmp.bed'
    sort_file   = 'process/sort.bed'
    merge_file  = 'process/merge.bed'
    merge_peaks = 'process/merged_enhancer_peaks_{}.srt.bed'.format(distance)

    if files_exist([merge_peaks]):
        return merge_peaks

    # Exclude peaks that overlap tss
    exclude_command = "bedtools intersect -a {} -b {} -v".format(peaks, tss)
    logger.info("Running exclude with %s", exclude_command)
    subprocess.check_call(exclude_command, stdout=open(tmp_file, 'w'), shell=True)

    # Sort and get unique rows
    sort_command = "bedtools sort -i {} | uniq > {}".format(tmp_file, sort_file)
    logger.info("Running exclude with %s", sort_command)
    subprocess.check_call(sort_command, shell=True)

    # Merge peaks
    merge_command = "bedtools merge -i {} -d {}".format(sort_file, distance)
    logger.info("Running merge with %s", merge_command)
    subprocess.check_call(merge_command, stdout=open(merge_file, 'w'), shell=True)

    # Sort merged file
    sort_command = "bedtools sort -i {} > {}".format(merge_file, merge_peaks)
    logger.info("Running sort with %s", sort_command)
    subprocess.check_call(sort_command, shell=True)

    os.remove(sort_file)
    os.remove(tmp_file)
    os.remove(merge_file)
    return merge_peaks

def bed_distribution(bed, fname, flabel):
    '''Explore the enhancer distribution on the merged peaks.'''

    enh = pd.read_csv(bed, sep='\t', header=None)
    distribution = 'results/{}_distribution_{}.png'.format(fname, flabel)

    # Plot enhancer disctribtion

    plt.hist(enh[2]-enh[1], bins=100)
    plt.xlabel('size (bp)')
    plt.ylabel('frequency')
    plt.savefig(distribution)
    plt.clf()

def enhancer_bounds(tss, 
                    region, 
                    sizes, 
                    distance, 
                    peaks, 
                    prange, 
                    logger):
    if len(prange) > 2:
        print("Please only use two numbers for Promoter range separated by a -",
              file = sys.stderr)
        exit(1)

    # download chromosome sizes
    if not os.path.exists(sizes):
        if not os.path.exists('process/%s.chrom.sizes' % sizes):
            os.system("wget hgdownload.cse.ucsc.edu/goldenPath/%s/"\
                            "bigZips/%s.chrom.sizes -P process/" % (sizes, sizes))
        sizes = 'process/%s.chrom.sizes' % sizes

    # Create a file handler
    handler = logging.FileHandler('log/enhancers_boundries.log')
    logger.addHandler(handler)

    # Run extend_tss
    extended_tss, promoter_file = extend_tss(tss, 
                                             sizes, 
                                             region, 
                                             prange, 
                                             logger)

    # Merge all peak files
    universe_peaks = merge_peaks(peaks[0], logger, None)
    peaks = peaks[1:]
    if len(peaks) > 1:
        for p in peaks:
            universe_peaks = merge_peaks(p, logger, universe_peaks)

    # Exlcude out tss and merge merge peaks
    merged = exclude_tss(universe_peaks, 
                         extended_tss, 
                         distance, 
                         logger)
    os.remove(universe_peaks)

    # Plot distribution
    bed_distribution(merged, 
                     'enhancer',
                     '_exclude_{}_within_{}'.format(region, distance))
    bed_distribution(promoter_file, 
                     'promoter', 
                     '_{}up_{}down'.format(prange[0], prange[1]))

    return merged, promoter_file

# --------------------------------------------------------------------------- #
# Methods for Merging Transcript Expression and getting a normalized log2     #
# based TPM matrix of expression                                              #
# --------------------------------------------------------------------------- #
def filter_txn_mx(tfmtx):
    corrmtx = np.square(tfmtx.corr())
    plot_corr_heatmap(corrmtx, "results/prefilter_tx_expression_log2_corr", 24)

    mask = corrmtx.columns[corrmtx.mean() > 0.5].to_list()
    tfmtx = tfmtx[mask]
    corrmtx = np.square(tfmtx.corr())
    plot_corr_heatmap(corrmtx, "results/postfilter_tx_expression_log2_corr", 24)

    return tfmtx

def get_filtered_transcripts(df, name):
    '''Filter for only ENS transcript.'''
    transcript_list     = list(df['gene_id'])
    selected_genes      = ['ENS' in c for c in transcript_list]
    selected_df         = pd.DataFrame(df[selected_genes][['gene_id', 'TPM']])
    selected_df.columns = ['gene_id', name]
    return selected_df

def merge_transcript(expression_files, logger):
    for fname in expression_files:
        if not os.path.exists(fname):
            print("Path %s does not exist!" % (fname),
                  file=sys.stderr)
            exit(1)

    # Create a file handler
    handler = logging.FileHandler('log/merge_transcript_expression.log')
    logger.addHandler(handler)

    # Get list of all expression matrixes
    list_files = []
    for files in expression_files:
        list_files.extend(glob.glob(files))

    # Loop through all files to add to dataframe
    expression_df = pd.DataFrame()
    for f in list_files:

        # Read in file name
        current_df = pd.read_csv(f, sep='\t')
        file_name  = os.path.basename(f)
        basename   = file_name.rpartition('.tsv')[0]
        logger.info('Processing: %s', f)

        # Filter genes
        selected_current_df = get_filtered_transcripts(current_df, basename)
        selected_current_df.set_index('gene_id', inplace=True)
        del(current_df)

        # Append to dataframe
        if expression_df.shape == (0, 0):
            expression_df = selected_current_df
        else:
            expression_df = expression_df.join(selected_current_df)

    # Convert TPM to log2(TPM)
    expression_nfile      = 'process/transcript_expression_log2TPM_matrix.txt'
    expression_df_log2TPM = np.log2(expression_df+1)

    expression_df_log2TPM = filter_txn_mx(expression_df_log2TPM)
    expression_df_log2TPM.to_csv(expression_nfile, sep='\t')

    plot_mtx(expression_nfile)
    return expression_nfile

# --------------------------------------------------------------------------- #
# Methods for Merging Enhancer Expression to get log2 TPM matrix              #
# --------------------------------------------------------------------------- #
def enhancer_coverage(alignment_tup):
    '''Get coverage for enhancers in sample.'''
    enhancer, alignment, logger = alignment_tup

    # Calculate coverage
    file_name        = os.path.basename(alignment)
    basename         = file_name.rpartition('.bam')[0]
    coverage_bed     = "process/{}_coverage.bed".format(basename)
    coverage_command = "bedtools coverage -sorted "\
                            "-a {} -b {} > {}".format(enhancer, 
                                                      alignment, 
                                                      coverage_bed)
    logger.info('Processing coverage using: %s', coverage_command)
    subprocess.check_call(coverage_command, shell=True)


def merge_enhancer(enhancers, target, alignment_files, logger):
    # Create a file handler
    handler = logging.FileHandler('log/merge_enhancer_expression.log')
    logger.addHandler(handler)

    # Get list of all alignment files
    list_files = []
    for files in alignment_files:
        list_files.extend(glob.glob(files))

    # Calculate coverage on all files
    all_files_tup = [(enhancers, x, logger) for x in list_files]
    cpus          = mp.cpu_count()
    process       = mp.Pool(cpus)

    process.map(enhancer_coverage, all_files_tup)
    process.close()
    process.join()

    # Loop through all files to add to dataframe
    enhancers_df          = pd.read_csv(enhancers, sep='\t', header=None)
    expression_df         = enhancers_df.iloc[:, 0:3]
    expression_df.columns = ['chr', 'start', 'end']
    del(enhancers_df)

    # List of coverage files
    list_of_coverage = glob.glob('{}/process/*_coverage.bed'.format(os.getcwd()))
    for f in list_of_coverage:
        # Read in file name
        file_name = os.path.basename(f)
        basename  = file_name.rpartition('_coverage.bed')[0]
        logger.info('Processing: %s', basename)

        # Calculate TPM
        counts_df       = pd.read_csv(f, sep='\t', header=None)
        enhancer_tpm    = (counts_df[3]/counts_df[5])\
                            *(1/sum(counts_df[3]/counts_df[5]))\
                            *1e6
        enhancer_signal = np.log2(enhancer_tpm+1)

        expression_df[basename] = enhancer_signal
        os.remove(f)

    # Convert TPM to log2(TPM)
    expression_df['name'] = expression_df['chr'] \
                                + ':' \
                                + expression_df['start'].map(str) \
                                + '-' \
                                + expression_df['end'].map(str)
    columns = ['chr', 'start', 'end']
    expression_df.drop(columns, inplace=True, axis=1)
    expression_df.set_index('name', inplace=True)
    expression_file="process/enhancer_{}_log2TPM_signal_matrix.txt".format(target)
    expression_df.to_csv(expression_file, sep='\t')

    bedenh = "process/enh.bed"
    os.system("awk '{print $1}' %s "\
                "| sed '1d' "\
                "| tr ':' '\\t' "\
                "| tr '-' '\\t' "\
                    "> %s" % (expression_file, bedenh))

    plot_mtx(expression_file)
    return expression_file, bedenh


# --------------------------------------------------------------------------- #
# Methods for Building Transcription Factor logodds matricies using FIMO      #
#   from meme suir of software                                                #
# --------------------------------------------------------------------------- #
def download_meme(database, species):
    if not os.path.exists("motif_databases"):
        cmdpipe = ["wget http://meme-suite.org/meme-software/"\
                        "Databases/motifs/motif_databases.12.18.tgz -P process/",
                   "tar xzf process/motif_databases.12.18.tgz"]

        for cmd in cmdpipe:
            try:
                # subprocess.check_call(shlex.split(cmd), shell=True)
                os.system(cmd)
            except:
                pass

    memels  = glob.glob("motif_databases/{}/*".format(database.upper()))
    species = '_'.join(species)
    memedir = "motif_databases/{}/{}.meme".format(database.upper(), species)
    if memedir not in memels:
        print("MEME database {} not found in {}".format(species, database))
        exit(1)
    
    return memedir

def download_genome(orgname, gvers, logger):
    path    = '/'.join(os.path.realpath(__file__).split('/')[:-2])
    cmd     = [path + '/trenco_tools/get_genome.py', '--refseq', '-s', orgname]
    genomes = subprocess.check_output(cmd + ['--get-versions'])
    genomes = genomes.decode("utf-8").replace('\n','\t').split('\t')
    genome  = ''
    print(genomes)
    for i in range(4,len(genomes),4):
        if gvers:
            if gvers in genomes[i]:
                genome = genomes[i]
        else:
            if genomes[i-1] == "latest":
                genome = genomes[i]

    assert genome
    if os.path.exists(genome + '_genomic.fna.gz'):
        subprocess.call(['gunzip', genome  + '_genomic.fna.gz'])
    elif os.path.exists(genome + '_genomic.fna'):
        genome = genome + '_genomic.fna'
    else:
        if gvers:
            cmd += ['-v', gvers]

        gout = subprocess.check_output(cmd + ['--use-ucsc'])
        if "Completed" not in str(gout):
            logger.warning('Error: Download output - {}'.format(gout))
            exit(1)
    
        genome = gout.decode("utf-8")\
                    .replace('\n','')\
                    .split(':')[-1]\
                    .replace('.gz','')
        
    return genome

def fimo_2_mtx(nfile, 
               name, 
               gff = False, 
               correction = None):
    assert os.path.exists(nfile)

    if '/' in name:
        name = name.split('/')[-1]

    tfmx = 'process/TFx{}_score_matrix.txt'.format(name)
    if os.path.exists(tfmx):
        return tfmx

    if gff: # NOT RECOMMENDED: the score is different
        melted_df = pd.read_csv(nfile, 
                                sep="\t", 
                                skiprows=1, 
                                header=None, 
                                engine='python')
        melted_df.columns = ['seqname', 
                             'source', 
                             'feature', 
                             'start', 
                             'end', 
                             'score', 
                             'strand', 
                             'frame', 
                             'attribute']

        # Function to parse the attribute column into multiple named 
        # columns of attributes VERY memory intensive
        attributes = pd.DataFrame.from_dict(
            [dict(
                a.split('=') for a in b if a
             ) for b in melted_df['attribute'].str.split(";")
            ])
        try:
            melted_df['MotifID'] = attributes['Alias']
        except:
            melted_df['MotifID'] = attributes['Name']

        melted_df = melted_df[['MotifID', 'seqname', 'score']]

    else:
        melted_df = pd.read_csv(nfile, sep="\t", skipfooter=4, engine='python')
        melted_df = melted_df[['motif_alt_id', 'sequence_name', 'score']]
        melted_df.columns = ['MotifID', 'seqname', 'score']

    #plot the histogram of the scores
    s = np.float64(melted_df[['score']])
    plt.hist(s, bins=100)
    plt.xlabel('FIMO score')
    plt.ylabel('frequancy')
    plt.savefig('results/fimo_score_freq_histi{}.png'.format(name))
    plt.clf()

    df_mtx = pd.pivot_table(melted_df, 
                            index='MotifID', 
                            columns='seqname', 
                            values='score')
    df_mtx.fillna(0, inplace=True)

    if correction:
        converter = load_convert_name(correction)
        df_mtx = df_mtx.rename(converter, axis=0)

    df_mtx.to_csv(tfmx, sep='\t')
    return tfmx

def run_fimo(odir, 
             genome, 
             motif_file, 
             bed, 
             logger, 
             correction = None, 
             wiggle = None, 
             use_priors = False):
    verscmd = ['bedtools', 
               '--version']
    bedvers = subprocess.check_output(verscmd, shell=True).decode('utf-8')
    bedvers = re.search(r'v([\d.]+)', bedvers).group(1).replace(".", "")

    cmdpipe = [['bedtools', 
                'getfasta', 
                '-fi', genome, 
                '-bed', bed, 
                '-fo', 
                '{}.fasta'.format(bed)]]
    if "promoter" in bed:
        if int(bedvers) < 2290:
            cmdpipe[-1].append("-name")
        else:
            cmdpipe[-1].append("-nameOnly")
    fullout = '{}/{}.motif'.format(odir, bed)

    if not os.path.exists(fullout):   
        if use_priors:
            cmdpipe += [] # TODO finish priors script
        else:
            cmdpipe += [['fimo', 
                         '-o', 
                         fullout, 
                         motif_file, 
                         '{}.fasta'.format(bed)]]

        files_needed = [genome, motif_file, bed]
        for f in files_needed:
            if not os.path.exists(f):
                logger.warning('Cannot find file: {}'.format(f))
                exit(1)

        for cmd in cmdpipe:
            print(" ".join(cmd))
            subprocess.check_call(cmd)

    fimo_fn = '{}/fimo.tsv'.format(fullout)
    tfmx = fimo_2_mtx(fimo_fn, bed.split('/')[-1], correction = correction)
    return tfmx


def build_tf_matrix(bed, 
                    meme_db, 
                    memedb, 
                    orgname, 
                    gvers, 
                    logger, 
                    vector = None):
    correction = check_if_commons(meme_db, False)
    meme_file  = download_meme(meme_db, memedb)
    genome     = download_genome(orgname, gvers, logger)
    fn_mx      = run_fimo('.', genome, meme_file, bed, logger, correction)

    logger.info("Using MEME: {}".format(meme_file))
    logger.info("Using Genome: {}".format(genome))
    logger.info("Built Motif file: {}".format(genome))
    
    return fn_mx


# --------------------------------------------------------------------------- #
# Main methods for getting Enhancer Gene Networks and building small TAD      #
#   matricies of enhancers x genes                                            #
# --------------------------------------------------------------------------- #
def init_e(_exp, _enh, _k):
    """To provide global variables shared across multiple processes."""
    global s_gene_tpm
    global s_enh_tpm
    global sample
    s_gene_tpm = _exp
    s_enh_tpm = _enh
    sample = _k

def get_distance_df(d):
    """To get distance between enh and gene for each TAD."""
    enh_path  = '{}/enh.bed'.format(d)
    gene_path = '{}/gene.bed'.format(d)
    if os.path.isfile(enh_path) and os.path.isfile(gene_path):
        enh_df  = pd.read_csv(enh_path, sep='\t', header=None)
        gene_df = pd.read_csv(gene_path, sep='\t', header=None)
        # generate enh identifiers
        enh_pos = []
        for r in enh_df.iterrows():
            pos = "{}:{}-{}".format(r[1][0], r[1][1], r[1][2])
            enh_pos.append(pos)
        # get gene names
        gene_names = list(gene_df[3])
        out_d = pd.DataFrame(0, index=enh_pos, columns=gene_names)
        # calculate enh gene distance
        for g in gene_df.iterrows():
            name = g[1][3]
            if g[1][5] == '-':
                g_start = g[1][2]
            else:
                g_start = g[1][1]
            for e in enh_df.iterrows():
                e_pos = "{}:{}-{}".format(e[1][0], e[1][1], e[1][2])
                e_mid = e[1][1]+((e[1][2]-e[1][1])/2)
                distance = np.abs(g_start - e_mid)
                out_d.loc[e_pos, name] = int(distance)

        out_d.to_csv('{}/enh_gene_distance.txt'.format(d), sep='\t')

        scaled_out_d = get_distance_based_weight(out_d)
        scaled_out_d.to_csv('{}/enh_gene_distance_scaled.txt'.format(d),
                            sep='\t')

def get_distance_based_weight(distance):
    """To scale the distance between enh and gene for each TAD."""
    # set weight 0 for distance greater than 500 kb
    # TODO clean this out of the code if we determine it to be unnecessary
    # distance[distance > 500000] = 1e10000
    dist_weight = 1/np.log2(distance)
    dist_list = np.sort(list(set(dist_weight.values.flatten())))

    if set(dist_weight.values.flatten()) == {0} or len(dist_list) == 1:
        scaled_dist_weight = dist_weight
    else:
        wmin = np.sort(list(set(dist_weight.values.flatten())))[1]
        wmax = np.sort(list(set(dist_weight.values.flatten())))[-1]
        scaled_dist_weight = (dist_weight-wmin)/(wmax-wmin)
        scaled_dist_weight[scaled_dist_weight < 0] = 0
    return(scaled_dist_weight)

def get_enh_gene_weights(d):
    """To get the enh gene weight based on enh activity and distance."""
    dist_path = '{}/enh_gene_distance_scaled.txt'.format(d)
    weight_path = '{}/{}_enh_gene_weight.txt'.format(d, sample)

    if os.path.exists(dist_path) and not os.path.exists(weight_path):
        dist_df = pd.read_csv(dist_path, sep='\t', index_col=0)
        tad_gene_tpm = s_gene_tpm.loc[
                            s_gene_tpm.index.intersection(dist_df.columns)
                        ].reindex(dist_df.columns)
        tad_enh_tpm = s_enh_tpm.loc[
                            s_enh_tpm.index.intersection(dist_df.index)
                        ].reindex(dist_df.index)

        try:
            enh_gene_weight = np.multiply(
                dist_df, np.sqrt(np.matmul(tad_enh_tpm, tad_gene_tpm.T)))
        except:
            enh_gene_weight = np.multiply(
                dist_df, np.sqrt(tad_enh_tpm.dot(tad_gene_tpm.T)))

        enh_gene_weight.to_csv(weight_path, sep='\t')
            
def plot_scatter(x_vals, y_vals, outdir, tad):
    #To scatter plot the scaled distance as function of enh gene distance.
    plt.scatter(x_vals, 
                y_vals, 
                edgecolors='#2374AB',
                facecolors='none',
                s=30,
                linewidths=1,
                alpha=0.4)
    plt.xscale('log')
    plt.vlines(500000, 0, 0.9, label='500kb')
    plt.xlabel('enh-gene distance (log scale)')
    plt.ylabel('weight')
    plt.legend()
    plt.savefig('{}/{}-scatter_enh_gene_scaled_distance.png'.format(outdir,tad))
    plt.clf()

def plot_heatmap(df, outdir, tad):
    """To plot heatmap of the scaled distance between enh and gene."""
    sb.set_context('paper', font_scale=1.8)
    g = sb.heatmap(df, yticklabels=False, cmap='Blues')
    fig = g.get_figure()
    fig.savefig('{}/{}-heatmap_enh_gene_scaled_distance.png'.format(outdir,tad))
    plt.clf()

def get_enh_networks(sample, 
                     gene, 
                     enh, 
                     sample_out, 
                     threads, 
                     enhBED, 
                     geneBED, 
                     tadBED, 
                     mask = []):
    os.makedirs(os.path.dirname('log/{}_log_enh_network.txt'.format(sample)),
                exist_ok=True)

    # get the enhancer and gene quantifications
    gene_tpm = pd.read_csv(gene, sep='\t', index_col=0)
    enh_tpm = pd.read_csv(enh, sep='\t', index_col=0)

    if mask:
        gene_mask = [s in mask for s in gene_tpm.columns]
        enh_mask = [s in mask for s in enh_tpm.columns]
    else:
        gene_mask = [sample in s for s in gene_tpm.columns]
        enh_mask = [sample in s for s in enh_tpm.columns]
    # Mean is done here
    s_gene_tpm         = pd.DataFrame(gene_tpm.loc[:,gene_mask].mean(axis=1)) 
    s_gene_tpm.columns = [sample]
    s_enh_tpm          = pd.DataFrame(enh_tpm.loc[:,enh_mask].mean(axis=1))
    s_enh_tpm.columns  = [sample]

    p = mp.Pool(threads, 
                initializer=init_e,
                initargs=(s_gene_tpm, s_enh_tpm, sample,))

    # Do the enh-gene distance calculation if not calculated
    if not os.path.exists('process/{}/TADs'.format(sample_out)) \
            or not os.path.exists('results/{}'.format(sample_out)):
        # ouput directories
        if not os.path.exists('process/{}'.format(sample_out)):
            os.system('mkdir process/{}'.format(sample_out))
        if not os.path.exists('results/networks'):
            os.system('mkdir results/networks')
        os.system('mkdir process/{}/TADs'.format(sample_out))
        os.system('mkdir results/{}'.format(sample_out))
    
        # sort tad bed file
        os.system('bedtools sort -i {} '\
                        '> process/{}/sort_enh.bed'.format(enhBED,
                                                           sample_out))
        os.system('bedtools sort -i {} '\
                        '> process/{}/sort_gene.bed'.format(geneBED,
                                                            sample_out))
        os.system('bedtools sort -i {} '\
                        '> process/{}/sort_tad.bed'.format(tadBED,
                                                           sample_out))

        ###############################
        # enh and gene overlapping tad
        print('Determining TAD overlapping enhancers and genes for %s' 
                % sample_out)
        os.system('bedtools intersect '\
                    '-a process/{a}/sort_enh.bed '\
                    '-b process/{a}/sort_tad.bed '\
                    '-wa -u -f 0.9 '\
                        '> process/{a}/sort_enh_in_tad.bed'.format(
                            a=sample_out
                        )
                 )
        os.system('bedtools intersect '\
                    '-a process/{a}/sort_gene.bed '\
                    '-b process/{a}/sort_tad.bed '\
                    '-wa -u -f 0.9 '\
                        '> process/{a}/sort_gene_in_tad.bed'.format(
                            a=sample_out
                        )
                  )
        os.system('bedtools intersect '\
                    '-a process/{a}/sort_enh_in_tad.bed '\
                    '-b process/{a}/sort_tad.bed '\
                    '-wao -f 0.9 '\
                        '> process/{a}/sort_enh_in_tad_detail.bed'.format(
                            a=sample_out
                        )
                  )
        os.system('bedtools intersect '\
                    '-a process/{a}/sort_gene_in_tad.bed '\
                    '-b process/{a}/sort_tad.bed '\
                    '-wao -f 0.9 '\
                        '> process/{a}/sort_gene_in_tad_detail.bed'.format(
                            a=sample_out
                        )
                 )
        ########################################
        # generate directory for each TAD region
        print('Generating directory for each TAD')
        tad_bed = pd.read_csv('process/{}/sort_tad.bed'.format(sample_out),
                              sep='\t', 
                              header=None)

        tad_dic = {}
        tad_key_dic = {}
        for indx, r in enumerate(tad_bed.iterrows()):
            pos = "{}:{}-{}".format(r[1][0], r[1][1], r[1][2])
            tad_dic[pos] = "TAD_{}".format(indx)
            tad_key_dic["TAD_{}".format(indx)] = pos

        for i in tad_dic.keys():
            os.system('mkdir process/{}/TADs/{}'.format(sample_out,
                                                        tad_dic[i]))

        print('Parsing enhancer and gene for each TAD')
        ###################
        # enh for each tad
        enh_tads = pd.read_csv(
            'process/{}/sort_enh_in_tad_detail.bed'.format(sample_out),
            sep='\t', 
            header=None
        )

        enh_tads_dic = {}
        for r in enh_tads.iterrows():
            tad_pos = "{}:{}-{}".format(r[1][3], r[1][4], r[1][5])
            enh_pos = (r[1][0], r[1][1], r[1][2])
            tad_key = tad_dic[tad_pos]
            if tad_key not in enh_tads_dic.keys():
                enh_tads_dic[tad_key] = [enh_pos]
            else:
                enh_tads_dic[tad_key].append(enh_pos)

        for k in enh_tads_dic.keys():
            enh_df = pd.DataFrame(enh_tads_dic[k])
            enh_df.to_csv('process/{}/TADs/{}/enh.bed'.format(sample_out, k),
                          sep='\t', 
                          index=None, 
                          header=None)

        ###################
        # gene for each tad
        gene_tads = pd.read_csv(
                    'process/{}/sort_gene_in_tad_detail.bed'.format(sample_out),
                    sep='\t', 
                    header=None)

        gene_tads_dic = {}
        for r in gene_tads.iterrows():
            tad_pos = "{}:{}-{}".format(r[1][6], r[1][7], r[1][8])
            gene_pos = (r[1][0], r[1][1], r[1][2], r[1][3], r[1][4], r[1][5])
            tad_key = tad_dic[tad_pos]
            if tad_key not in gene_tads_dic.keys():
                gene_tads_dic[tad_key] = [gene_pos]
            else:
                gene_tads_dic[tad_key].append(gene_pos)

        for k in gene_tads_dic.keys():
            gene_df = pd.DataFrame(gene_tads_dic[k])
            gene_df.to_csv('process/{}/TADs/{}/gene.bed'.format(sample_out, k),
                           sep='\t', 
                           index=None, 
                           header=None)

        print('Calculating enhancer and gene distance for each TAD')
        # calculate enh-gene distance for each TAD
        tad_dirs = glob.glob('process/{}/TADs/TAD_*'.format(sample_out))
        out      = p.map(get_distance_df, tad_dirs)

        # plot enh gene distance weights
        for tad in tad_dirs:
            dist_path   = '{}/enh_gene_distance.txt'.format(tad)
            s_dist_path = '{}/enh_gene_distance_scaled.txt'.format(tad)

            if not os.path.exists(dist_path) or not os.path.exists(s_dist_path):
                continue

            dist        = pd.read_csv(dist_path, sep='\t', index_col=0)
            scaled_dist = pd.read_csv(s_dist_path, sep='\t', index_col=0)
            tadfn       = tad.split('/')[-1]
            plot_scatter(dist.values.flatten(), 
                         scaled_dist.values.flatten(),
                        'results/{}'.format(sample_out), 
                        tadfn)
            plot_heatmap(scaled_dist, 
                         'results/{}'.format(sample_out), 
                         tadfn)

    print('Calculating enhancer gene weights for each TAD')
    # calculate enh-gene weights for each TAD
    tad_dirs = glob.glob('process/{}/TADs/TAD_*'.format(sample_out))

    if s_gene_tpm.shape[1] == 1 and s_enh_tpm.shape[1] == 1:
        p.map(get_enh_gene_weights, tad_dirs)
        p.close()
    else:
        print('There are multiple samples with same name. Please supply' \
               'dataframe containing columns with unique names',
               file=sys.stderr)
        exit(1)

    gene_vec_name = 'process/{}_gene_vector.txt'.format(sample)
    s_gene_tpm.to_csv(gene_vec_name, sep='\t')

    return gene_vec_name

# --------------------------------------------------------------------------- #
# Main methods for generating full Tissue Specific Gene Network               #
# --------------------------------------------------------------------------- #
def init(_s, _m):
    global sample
    global mx
    sample = _s
    mx = _m

def get_expectation(d):
    """To get expectation of the networks for each TAD."""
    expectation_file = '{}/expectation.txt'.format(d)
    if not os.path.exists(expectation_file):
        tad_files = glob.glob('{}/*_enh_gene_weight.txt'.format(d))
        if not tad_files:
            print("No gene_weight file in {}".format(d))
            return()

        # calculate the expectation
        for indx, k in enumerate(tad_files):
            if indx == 0:
                sum_enh_gene = pd.read_csv(k, sep='\t', index_col=0)
            else:
                other_enh_gene = pd.read_csv(k, sep='\t', index_col=0)
                sum_enh_gene += other_enh_gene
        expectation = sum_enh_gene/len(tad_files)
        expectation.to_csv('{}/expectation.txt'.format(d), sep='\t')

def get_transformed_networks(d):
    """To generate tissue specific network for each TAD."""
    try:
        expectation = pd.read_csv('{}/expectation.txt'.format(d), sep='\t',
                                  index_col=0)
    except Exception:
        return() 
    network = pd.read_csv('{}/{}_enh_gene_weight.txt'.format(d, sample),
                          sep='\t', index_col=0)
    w = np.arctan(network - expectation) # look into how arctan will work with this
    w.to_csv('{}/{}_enh_gene_weight_minus_exp_tanh.txt'.format(d, sample),
             sep='\t')

def get_network(d):
    """To get the network for each TAD."""
    enh_path     = '{}/enh.bed'.format(d)
    # Only using weight here rather than enh_gene_weight_minus_exp_tanh.txt
    weight_path  = '{}/{}_enh_gene_weight.txt'.format(d, sample) 
    tf_gene_path = '{}/tf_gene_wgt_matrix.txt'.format(d)
    if os.path.exists(enh_path) \
            and os.path.exists(weight_path) \
            and not os.path.exists(tf_gene_path):
        weight = pd.read_csv(weight_path, sep='\t',index_col=0)
        enh_list = list(weight.index)

        enh_mask = [x in enh_list for x in list(mx.columns)]
        mx_sub = pd.DataFrame(mx.loc[:,enh_mask])
        mx_sub.to_csv('{}/tf_enh_matrix.txt'.format(d), sep='\t')

        try:
            weight_mx = np.matmul(mx_sub, weight)
        except:
            weight_mx = mx_sub.dot(weight)

        weight_mx.to_csv(tf_gene_path, sep='\t')

        return weight_mx

    elif os.path.exists(tf_gene_path):
        weight_mx = pd.read_csv(tf_gene_path, sep='\t',index_col=0)
        return weight_mx

    else:
        return None

def tis_enh_gene_net(sample, 
                     outdir, 
                     matrix1, 
                     matrix2, 
                     gene_vect, 
                     threads):
    if not os.path.exists(matrix1) \
            or not os.path.exists(matrix2) \
            or not os.path.exists(gene_vect):
        print("Matrix file or Gene Vector file cannot be found",
              file=sys.stderr)
        exit(1)

    # Record matrix and matrix distribution
    def record_mx(mtx, name):
        mtx.to_csv('process/{}_{}.txt'.format(sample, name), sep='\t')

        # plot histogram distribution
        flat_mtx = mtx.values.flatten()
        pt       = sb.distplot(flat_mtx, hist_kws={'log':True})
        fig      = pt.get_figure()
        fig.savefig("results/%s/%s.png" % (sample, name))
        plt.clf()

    # get the TF-by-Enh matrix and TF-by-Promoter
    tf_enh_mx = pd.read_csv(matrix1, sep='\t', index_col=0)
    tf_pro_mx = pd.read_csv(matrix2, sep='\t', index_col=0)
    g_vect    = pd.read_csv(gene_vect, sep='\t', index_col=0)

    # TFscale = 40 # add a scale to the TF matricies that compresses 
    # the data into a 0-1 range
    tf_enh_mx = tf_enh_mx / tf_enh_mx.values.max()
    tf_pro_mx = tf_pro_mx / tf_pro_mx.values.max()

    # This filtering will set the g_vector to the selected regions 
    # (ie protein coding)
    # we can consider expanding the selection to and not allowing selection 
    # of regions
    inter     = g_vect.index[[s in tf_pro_mx.columns for s in g_vect.index]]
    g_vect    = pd.DataFrame(g_vect.loc[inter,:])
    tf_pro_mx = pd.DataFrame(tf_pro_mx.loc[:,inter])

    # Tissue specific enh gene networks by subtracting the expectation
    print('Calculate weights of Enh-Gene by subtracting expectation')
    tad_dirs = glob.glob('process/{}/TADs/TAD_*'.format(outdir))
    p = mp.Pool(threads,
                initializer=init,
                initargs=(sample, tf_enh_mx))

    dfs = p.map(get_network, tad_dirs)

    print('Number of individual TAD network', len(dfs))
    p.close()
    p.join()

    print('Gathering all weights from all TADs')  
    # Get the main matricies that we need  
    tf_gene_mx_pro = tf_pro_mx.multiply(g_vect.T.iloc[0,:], axis=1)
    record_mx(tf_gene_mx_pro, "tf_gene_promoter_matrix")
    
    tf_gene_mx_enh = pd.concat(dfs, sort = False, axis = 1)
    tf_gene_mx_enh = tf_gene_mx_enh.loc[:,inter].fillna(0)
    record_mx(tf_gene_mx_enh, "tf_gene_enhancer_matrix")

    # transform matricies with arctan and median
    tf_gene_mx_pro = np.arctan(tf_gene_mx_pro)
    record_mx(tf_gene_mx_pro, "tf_gene_promoter_matrix_arctan") 

    tf_gene_mx_enh = np.arctan(tf_gene_mx_enh)
    record_mx(tf_gene_mx_enh, "tf_gene_enhancer_matrix_arctan")

    # load name conversion tables
    with open("process/converter.json", "r") as ifi:
        conversion = json.load(ifi)

    # Add matricies and build full network
    g_vect = g_vect.rename(conversion, axis=0)
    g_vect = g_vect.groupby(g_vect.index).mean()

    all_network  = tf_gene_mx_pro.add(tf_gene_mx_enh)/np.pi
    all_network  = all_network.rename(conversion, axis = 1)
    all_network  = all_network.groupby(all_network.index).mean()
    all_network  = all_network.groupby(all_network.columns, axis = 1).mean()
    intersection = all_network.index & all_network.columns
    all_network  = all_network.loc[intersection,:]

    tf_wgt = g_vect.loc[intersection,:]
    tf_wgt = tf_wgt / tf_wgt.max()
    record_mx(tf_wgt, "tfgene_weights_intermentiate_vector")

    all_network = all_network.mul(tf_wgt.iloc[:,0], axis=0)
    record_mx(all_network, "tf_gene_full_matrix")

    # Convert to vector based graph
    edges = all_network.stack().reset_index()
    edges.columns = ["TF", "Gene", "Weight"]
    edges = edges.sort_values("Weight", ascending=False)
    edges = edges.loc[(edges["Weight"] != 0.0),:]

    edges.to_csv('results/networks/{}_edges.csv'.format(sample), 
                 sep='\t', 
                 index = False)
    g_vect.to_csv('results/networks/{}_nodes.csv'.format(sample), sep='\t')

    print(all_network.shape)
    print('======================================')
