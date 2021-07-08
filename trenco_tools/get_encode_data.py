#!/usr/bin/env python3

'''Download data (ChIP-seq,RNA-seq) from ENCODE Portal using metadata.tsv.'''

import argparse
import pandas as pd
import os
import io
import requests
import multiprocessing as mp
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from hashlib import md5
import gzip
import logging

EPILOG = '''
For more details:
        %(prog)s --help
'''

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = True

_ASSEMBLY_MAPPER = {
    'hg38': 0,
    'GRCh38': 0,
    'hg19': 0,
    'mm10': 1,
    'mm9': 1,
    'dm6': 3,
    'dm3': 3,
    'ce11': 2,
    'ce10': 2,
    'ce6': 2,
}

_ASSEMBLY_KEY = ['Homo_sapiens', 'Mus_musculus', 'Caenorhabditis_elegans', 'Drosophila_melanogaster']

def get_args():
    '''Define arguments.'''

    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-a', '--assembly',
                        help="The data assembly to download. (e.g. mm10, hg19, GRCh38, etc.",
                        required=True)

    parser.add_argument('-o', '--out',
                        help="The root out directory to download files.",
                        required=True)

    parser.add_argument('-b', '--biosample',
                        help="The biosamples to filter for. (e.g. primary cell)",
                        nargs='+',
                        required=True)

    parser.add_argument('-t', '--age',
                        help="The ages to filter out. (e.g. 24w)",
                        nargs='+',
                        required=False)

    parser.add_argument('-s','--stage',
                        help="The life stages to filter out. (e.g. adult)",
                        nargs='+',
                        required=False)

    parser.add_argument('-p', '--proxy',
                        help="True/False if proxy required (UTSW).",
                        default=False,
                        action='store_true')

    args = parser.parse_args()
    return args


def get_metadata(assembly, organism, out_dir, proxy=False):
    '''Get the metadata.tsv file for all released Experiments.'''

    genus, species = organism.split('_')
    # Define metadata url
    metadata_url = "https://www.encodeproject.org/metadata/type=Experiment&status=released&assembly=%s/metadata.tsv" % (assembly)
    auxilary_url = 'https://www.encodeproject.org/report.tsv?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=%s+%s&field=accession&field=replicates.library.biosample.age&field=replicates.library.biosample.age_units&field=replicates.library.biosample.life_stage&assembly=%s' % (genus, species, assembly)

    # Proxy settings for UTSW
    if proxy:
        http_proxy = "http:proxy.swmed.edu:3128"
        https_proxy = "http:proxy.swmed.edu:3128"
        proxy_dict = {
                      "http": http_proxy,
                      "https": https_proxy
                    }
        metadata_file = requests.get(metadata_url, proxies=proxy_dict).content
        auxilary_file = requests.get(auxilary_url, proxies=proxy_dict).content
    else:
        metadata_file = requests.get(metadata_url).content
        auxilary_file = requests.get(auxilary_url).content

    # Date of file download
    date_download = datetime.datetime.now().strftime("%Y%m%d")
    download_file = "{}/metadata_{}_{}.tsv".format(out_dir,
                                                    assembly,
                                                    date_download)

    metadata_df = pd.read_csv(io.StringIO(metadata_file.decode('utf-8')),
                                sep='\t', low_memory=False)

    auxilary_df = pd.read_csv(io.StringIO(auxilary_file.decode('utf-8')),
                              sep='\t', 
                              low_memory=False, 
                              header = [1])

    merged_df = pd.merge(metadata_df, auxilary_df, left_on = ['Experiment accession'], right_on = ['Accession'])

    # Save file to directory
    merged_df.to_csv(download_file, sep='\t', encoding='utf-8', index=False)

    return merged_df


def convert_age(age, units, stage, organism):
    '''Return age in same format for converstion.'''

    converted_age = ['','']
    if age not in ['unknown']:
        if units == 'week':
            converted_age = ['{}w'.format(age), stage]
        elif units == 'day':
            converted_age = ['{}d'.format(age), stage]
        elif units == 'year':
            converted_age = ['{}y'.format(age), stage]
    else:
        converted_age = [age, stage]

    return converted_age


def get_df(df, assay_target, biosample_types, filter_age, output_type, file_type, assembly, organism):
    '''Return files matching filters.'''

    meta_tups = []
    df_filter = df[(df['Output type'] == output_type) & (df['File format']== file_type) & (df['Assembly']== assembly)]
    
    #os.mkdir('tmp/')
    #download_file = 'tmp/df_filter.csv'
    #df_filter.to_csv(download_file, sep='\t', encoding='utf-8', index=False)

    # Interate through rows to exlcude out data
    for r in df_filter.iterrows():

        # Define biosample varibles
        biosample_type = r[1]['Biosample type']
        biosample_term = '-'.join(r[1]['Biosample term name'].split(' '))
        biosample_age = r[1]['Age']
        biosample_age_units = r[1]['Age units']
        biosample_stage = r[1]['Life stage']

        # Define Experiment varibles
        bio_rep = r[1]['Biological replicate(s)']
        exp_acc = r[1]['Experiment accession']
        file_acc = r[1]['File accession']
        url = r[1]['File download URL']
        output = r[1]['Output type']
        status = r[1]['File Status']
        file_format = r[1]['File format']
        audit = r[1]['Audit ERROR']
        assay = r[1]['Assay']
        md5sum = r[1]['md5sum']

        # Check biosample type
        if (biosample_type not in biosample_types) and (biosample_term not in biosample_types):
            continue

        if (biosample_term in biosample_types):
            term = True
        else:
            assert (biosample_type in biosample_types)
            term = False

        # Set age
        age = convert_age(biosample_age, biosample_age_units, biosample_stage, organism)
        
        # Filter ages
        if filter_age[0] and filter_age[1]:
            if age[0] in filter_age[0] and age[1] in filter_age[1]:
                continue
        elif filter_age[0] or filter_age[1]:
            if age[0] in filter_age[0] or age[1] in filter_age[1]:
                continue

        # Filter out low read depth
        if audit in ["extremely low read depth"]:
            continue

        # Filter for released files
        if status != 'released':
            continue

        # Define Experiment target
        if assay_target == 'RNA-seq':
            target = assay
        elif assay == 'ChIP-seq':
            try:
                target = r[1]['Experiment target'].split('-')[0]
            except:
                continue
                pass
        else:
            continue

        # Make Tuple
        if target == assay_target:
            bio_rep = bio_rep.replace(', ', ',', 1)
            experiment = '{}_{}_rep{}_{}_{}_{}'.format(
                                                    biosample_term,
                                                    '-'.join(age),
                                                    bio_rep,
                                                    target,
                                                    file_acc,
                                                    exp_acc)

            term_age = '{}_{}'.format(biosample_type if term else biosample_term, age)
            output_type = output_type.replace(' ', '_')
            meta_tups.append(
                (target,
                url,
                file_acc,
                bio_rep,
                status,
                biosample_term if term else biosample_type,
                experiment,
                term_age,
                output_type,
                md5sum))

    # Return Dataframe of files matching criteria
    meta_columns = ['Target', 'URL', 'File Accession', 'Rep', 'Status', 'Type', 'Experiment', 'Age_tissue', 'Output Type', 'MD5Sum']
    
    filtered_df = pd.DataFrame(meta_tups, columns=meta_columns)
    if filtered_df.empty:
        print("No data found with given parameters:")
        print(assay_target, biosample_types, filter_age, output_type, file_type, assembly, organism)
        exit(1)
    return filtered_df

def overlap_experiments(peaks, alignments, quantifications):
    '''Return filtered dataframes with Experiemnts that have ChIP-seq and RNA-seq.'''

    # Define Unqiue sets
    unique_align = set(alignments['Age_tissue'])
    unique_peaks = set(peaks['Age_tissue'])
    unique_quant = set(quantifications['Age_tissue'])

    # Intersect ChIP-seq Alignments and Peaks
    unique_chip = set(unique_peaks & unique_align)

    # Intersect RNA-seq and ChIP-seq
    unique_experiments = set(unique_chip & unique_quant)

    # Filter for rows only matching Unique Experiments
    alignments_filter = alignments[alignments['Age_tissue'].isin(unique_experiments)]
    peaks_filter = peaks[peaks['Age_tissue'].isin(unique_experiments)]
    quantifications_filter = quantifications[quantifications['Age_tissue'].isin(unique_experiments)]

    return peaks_filter, alignments_filter, quantifications_filter


def avail_data(experiments, out_dir, biosample_type):
    '''Show plots of available data based on filters.'''

    cell_types = set([i.split('_', 1)[0] for i in experiments])
    age = sorted(set([i.split('_', 1)[1] for i in experiments]))
    experiment_df = pd.DataFrame(0, index=cell_types, columns=age)

    for i in experiments:
        cell_type = i.split('_', 1)[0]
        age = i.split('_', 1)[1]
        experiment_df.loc[cell_type, age] = 1

    file_name = '/{}_experiment_heatmap.png'.format(biosample_type)
    sns.set_context('paper', font_scale=1.8)
    hmap = sns.heatmap(experiment_df, vmin=0, vmax=1, cmap="Blues")
    plt.yticks(rotation=0)
    plt.xticks(rotation=45)
    plt.savefig(out_dir + file_name)
    plt.clf()

    return experiment_df

def getmd5(md5, filename):
    m = md5()
    p = gzip.open(filename, 'r+')
    m.update(p.read())
    p.close()
    checksum = m.hexdigest()
    return  checksum

def get_file(file_tup):
    '''Download files and rename appropriatley.'''

    mark, url, name, output_type, md5sum, out_dir = file_tup

    file_name = os.path.basename(url)
    extension = '.'.join(file_name.split(".")[1:])

    if mark != 'RNA-seq':
        assay_dir = '{}/{}/{}'.format(out_dir, 'ChIP-seq', mark)
    else:
        assay_dir = '{}/{}'.format(out_dir, mark)

    file_type_dir = '{}/{}'.format(assay_dir, output_type)

    lock.acquire()
    if not os.path.exists(assay_dir):
        os.makedirs(assay_dir)

    if not os.path.exists(file_type_dir):
        os.makedirs(file_type_dir)
    lock.release()

    # Download file and rename
    
    os.system('wget -P {} {}'.format(file_type_dir, url))
    os.system('mv {}/{} {}/{}.{}'.format(file_type_dir, file_name, file_type_dir, name, extension))

    # Check md5sum of content TODO: Fix md5sum check
    #ilename = '{}/{}.{}'.format(file_type_dir, name, extension)
    #download_md5sum = getmd5(md5, filename)
    #if download_md5sum != md5sum:
    #    logger.error("md5sum doesn't match for file: %s from %s" % (filename, url))


def main():
    args = get_args()
    assembly = args.assembly
    out =  args.out
    proxy = args.proxy
    biosample_types = args.biosample
    age = [[],[]]
    if args.age:
        if len(args.age) == 1 and args.age[0][0] in ['<','>']:
            if '<' == args.age[0][0]:
                suffix = args.age[0][-1]
                age[0] = [str(x)+suffix for x in range(0,int(args.age[1:-1]))]
            else:
                assert '>' == args.age[0][0]
                suffix = args.age[0][-1]
                age[0] = [str(x)+suffix for x in range(int(args.age[1:-1]),365)]
        else:
            age[0] = args.age
    if args.stage:
        age[1] = args.stage

    if age[1]:
        for i in range(len(age[1])):
            if age[1][i] in ['P','p']:
                age[1][i] = 'Postnatal'
            elif age[1][i] == 'e':
                age[1][i] = 'embryonic'

    # Check if assembly is valid for encode
    if assembly not in _ASSEMBLY_MAPPER.keys():
        raise Exception('Incorrect assembly')
    else:
        organism = _ASSEMBLY_KEY[_ASSEMBLY_MAPPER[assembly]]

    # Make out directory if not exists
    out_dir = '{}/{}/{}'.format(out, organism, assembly)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Create a file handler
    handler = logging.FileHandler('{}/download.log'.format(out_dir))
    logger.addHandler(handler)

    # Get metadata file
    metadata_df = get_metadata(assembly, organism, out_dir, proxy)
    
    # Filter for ChIP-seq H3K27ac Alignments and Replicated Peaks
    chip_align_df = get_df(metadata_df, 'H3K27ac', biosample_types, age, 'alignments', 'bam', assembly, organism)
    chip_peaks_df = get_df(metadata_df, 'H3K27ac', biosample_types, age, 'replicated peaks', 'bed narrowPeak', assembly, organism)

    # Filter for RNA-seq Transcript quantifications
    rna_quant_df = get_df(metadata_df, 'RNA-seq', biosample_types, age, 'gene quantifications', 'tsv', assembly, organism)

    # Filter for overlap experiments
    chip_peaks_filter_df, chip_align_filter_df, rna_quant_filter_df = \
                overlap_experiments(chip_peaks_df, chip_align_df, rna_quant_df)

    # Make heatmaps of avaialable data by Age and Tissue
    for t in biosample_types:
        experiments_df = chip_peaks_filter_df[chip_peaks_filter_df['Type'] == t]
        avail_data(set(experiments_df['Age_tissue']), out_dir, t)

    # Merge all dataframes
    all_files = pd.concat(
                    [chip_peaks_filter_df,
                    chip_align_filter_df,
                    rna_quant_filter_df])

    all_files['Out Dir'] = out_dir
    all_files_filter = all_files[['Target', 'URL', 'Experiment', 'Output Type', 'MD5Sum','Out Dir']]

    files = "{}/downloaded_files.tsv".format(out_dir)

    all_files_filter.to_csv(files, sep='\t', encoding='utf-8', index=False)

    all_files_tup = [tuple(x) for x in all_files_filter.to_records(index=False)]

    # Download all files
    def init(l):
        global lock
        lock = l
    cpus = mp.cpu_count()
    l = mp.Lock()
    process = mp.Pool(cpus, initializer=init, initargs=(l,))
    process.map(get_file, all_files_tup)
    process.close()
    process.join()

if __name__ == '__main__':
    main()
