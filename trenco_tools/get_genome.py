#!/usr/bin/env python3

import sys, os, subprocess, re
from argparse import ArgumentParser

_OMIT_COLUMNS = ['biosample', 
                 'wgs_master', 
                 'refseq_category', 
                 'infraspecific_name', 
                 'isolate', 
                 'relation_to_type_material',
                 'organism_name']

_DB_SOURCE_KEY = {'GCF' : 'RefSeq',
                  'GCA' : 'GenBank'}

# Read genome fasta files into memory
def read_genome(genome_file):
    chr_dic, chr_names, chr_full_names = {}, [], []

    seqs = open(genome_file, 'r').read()
    seqs = seqs.strip('\n').split('>')[1:]

    chr_name, chr_full_name, sequence = "", "", ""
    while seqs:
        ix = seqs[0].find('\n')

        chr_full_name, sequence = seqs[0][:ix], seqs[0][ix:].replace('\n','')
        chr_name = chr_full_name.split()[0]

        chr_dic[chr_name] = sequence
        chr_names.append(chr_name)
        chr_full_names.append(chr_full_name)

        seqs.pop(0)
    
    return chr_dic, chr_names, chr_full_names

# Write Sequence Dictionary to Fasta
def write_fasta(file_path, nfile, seqdict, split = True):
    if not os.path.exists(file_path) and file_path:
        os.makedirs(file_path)
    
    if file_path and not file_path.endswith('/'):
        file_path += '/'

    filename = file_path + nfile
    ofile = open(filename, 'w')
    
    for seq_name, seq in seqdict.items():
        ofile.write('>' + seq_name + '\n')

        if split:
            for subseq in (seq[i:i+80] for i in range(0, len(seq), 80)):
                ofile.write(subseq + '\n')
        else:
            ofile.write(seq + '\n')

    ofile.close()
    return filename

# Check if the database exists and if not build it from NCBI data
# Everything is positionally indexed
def download_index(dbname):
    if os.path.exists(dbname):
        return 0
    
    os.mkdir('tmp_idx')
    cmds = ['wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -P tmp_idx/',
            'wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt -P tmp_idx/',
            'wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -P tmp_idx/',
            'wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt -P tmp_idx/']
    
    files = []
    for cmd in cmds:
        cmd_ls = cmd.split()
        nfile = cmd_ls[1].split('/')[-1]
        files.append('tmp_idx/' + nfile)

        subprocess.run(cmd_ls)

    dataset = {}
    data_idx, col_count = {'bitmarker' : []}, 1
    db_colnames = ['bitmarker']
    for file_ in files:
        file_index, file_col = {}, {}
        with open(file_, "r") as ifi:
            for line in ifi:
                if "README" in line:
                    continue
                
                line = line.strip().split('\t')
                if line[0].startswith('#'):
                    line[0] = line[0].split()[-1]
                    for i in range(len(line)):
                        entry = line[i]
                        file_col[i] = entry
                        file_index[entry] = i
                        if entry in _OMIT_COLUMNS:
                            continue

                        if entry not in data_idx:
                            data_idx[entry] = col_count
                            db_colnames.append(entry)
                            col_count += 1
                    continue
                        
                assert file_index

                org_name = line[file_index['organism_name']]
                if org_name not in dataset:
                    dataset[org_name] = [[] for _ in range(col_count)]
                elif len(dataset[org_name]) < col_count:
                    dataset[org_name] = dataset[org_name] + [[] for _ in range(col_count - len(dataset[org_name]))]

                if line[file_index['ftp_path']] == 'na':
                    continue

                for i in range(len(line)):
                    fcol_nam = file_col[i]
                    if fcol_nam in _OMIT_COLUMNS:
                        continue
                
                    if i == 0:
                        if _DB_SOURCE_KEY[line[i][0:3]] == 'RefSeq':
                            baseid = '1'
                        else:
                            baseid = '0'
                
                    dbidx = data_idx[fcol_nam]
                    dataset[org_name][dbidx].append(line[i])

                dataset[org_name][0].append(baseid)

    with open(dbname, 'w') as ofo:
        ofo.write('# organism\t' + '\t'.join(db_colnames) + '\n')
        for org, line in dataset.items():

            newline = '\t'.join([','.join(sublist) for sublist in line])
            ofo.write(org + '\t' + newline + '\n')

    os.system('rm -r tmp_idx/')

# Load NCBI database of genomes that contains the FTP address and all information needed for script
def load_database(dbname):
    if not os.path.exists(dbname):
        print("Database not found. There may be a problem in the build")
        exit(1)

    db = {}
    with open(dbname, 'r') as ifi:
        for line in ifi:
            line = line.strip().split('\t')
            if line[0][0] == '#':
                dbcols = line
                continue
            
            assert dbcols

            org_name = line[0]
            db[org_name] = {}
            for i in range(1, len(line)):
                db[org_name][dbcols[i]] = line[i].split(',')

    return db

# NCBI genomes have chromosomes named as the contig ID they derive from.
# This script will convert those IDs to names
def convert_names(fasta, index, use_ucsc):
    whole_genome, chr_names, full_chr_names = read_genome(fasta)

    # read index of names
    ref2gen = {}
    gname_ix = {}
    with open(index, 'r') as ifi:
        for line in ifi:
            if line.startswith('#'):
                continue

            seqname, seqroll, asmol, astype, gen, rel, ref, asunit, seqlen, uscs = line.strip().split('\t')
            ref2gen[ref] = gen
            gname_ix[gen] = {'chr' : seqname,
                             'UCSC': uscs}

    for i in range(len(chr_names)):
        chr_name = chr_names[i]
        if chr_name in ref2gen:
            chr_name = ref2gen[chr_name]

        if use_ucsc:
            chr_name = gname_ix[chr_name]['UCSC']
        else:
            chr_name = gname_ix[chr_name]['chr']

        chr_name = chr_name + ' ' + full_chr_names[i]

        whole_genome[chr_name] = whole_genome.pop(chr_names[i])

    out_dir = '/'.join(fasta.split('/')[:-1])
    file_name = fasta.split('/')[-1]

    os.remove(fasta)
    write_fasta(out_dir, file_name, whole_genome)

if __name__ == '__main__':
    parser = ArgumentParser(
    description = "Get Genome by Species Name")

    parser.add_argument('-s', '--species',
                        dest='species',
                        nargs='+',
                        help = "Scientific name of organism (Can be left empty if using --list-all")
    parser.add_argument('--refseq',
                        dest='refseq',
                        action='store_true',
                        help = "Get RefSeq release of Genome")
    parser.add_argument('--genbank',
                        dest='genbank',
                        action='store_true',
                        help = "Get GenBank release of Genome")
    parser.add_argument('--get-annotations',
                        dest='get_ann',
                        action='store_true',
                        help = "Get the annotation file instead of the genome")
    parser.add_argument('-v', '--version',
                        dest='gvers',
                        type=str,
                        help = "Version of Genome to download. Leave blank for newest")
    parser.add_argument('-o', '--out-dir',
                        dest='out_dir',
                        type=str,
                        default="./",
                        help = "Directory to download Genome")
    parser.add_argument('--info',
                        dest='info',
                        action='store_true',
                        help = "Get information about the genome and version selected")
    parser.add_argument('--no-name-convert',
                        dest='dont_convert',
                        action='store_true',
                        help = "Convert Default NCBI name to standard name (ex. chr1, etc) ")
    parser.add_argument('--use-ucsc',
                        dest='use_ucsc',
                        action='store_true',
                        help = "Use UCSC naming if converting NCBI name (UCSC: chr1 vs 1)")
    parser.add_argument("--list-all",
                        dest='getall',
                        action='store_true',
                        help="List all species in database")
    parser.add_argument("--get-versions",
                        dest='getver',
                        action='store_true',
                        help="Get common availible versions for the species")

    args = parser.parse_args()

    if args.species:
        args.species = ' '.join(args.species)
    real_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    dbname = real_path + '/' + 'genome_database.txt'
    if args.out_dir[-1] != '/':
        args.out_dir += '/'

    download_index(dbname)

    genomes = load_database(dbname)
    species_list = list(genomes.keys())

    if args.getall:
        print(',\n'.join(species_list))
        exit(0)

    if args.species not in species_list:
        print("Can't locate species %s in database:" % args.species)
        possible = []
        for spec in species_list:
            argspec_ = args.species.split()
            for as_ in argspec_:
                if as_ in spec:
                    possible.append(spec)

        if not possible:
            print("\t0 possible matches")
        else:
            print("\tDid you perhaps mean:")
            for poss in possible:
                print("\t\t" + poss)

        exit(1)

    bitline = [True if x == '1' else False for x in genomes[args.species]['bitmarker']]
    asm_vers = genomes[args.species]['asm_name']
    if args.getver:
        new_asm = []

        for i in range(len(asm_vers)):
            status = genomes[args.species]['version_status'][i]
            nfile = genomes[args.species]['ftp_path'][i].split('/')[-1]

            if args.refseq and bitline[i]:
                new_asm.append('{}\trefseq\t{}\t{}'.format(asm_vers[i], status, nfile))
            if args.genbank and not bitline[i]:
                new_asm.append('{}\tgenbank\t{}\t{}'.format(asm_vers[i], status, nfile))

        new_asm = '\n'.join(set(sorted(new_asm)))
        print("{}:\n{}".format(args.species, new_asm))
        exit(0)

    if args.gvers and args.gvers not in asm_vers:
        print("{} genome version {} not found.".format(args.species, args.gvers))
        exit(1)

    indecies, base, websites = [], [], []
    for i in range(len(asm_vers)):
        if not args.gvers:
            if genomes[args.species]['version_status'][i] != 'latest':
                continue
            else:
                args.gvers = asm_vers[i]

        if args.gvers == asm_vers[i]:
            if args.refseq and bitline[i]:
                indecies.append(i)
                base.append('RefSeq')
                websites.append(genomes[args.species]['ftp_path'][i])
            if args.genbank and not bitline[i]:
                indecies.append(i)
                base.append('GenBank')
                websites.append(genomes[args.species]['ftp_path'][i])

    if not args.gvers:
        print('No annotated latest release found for %s!' % args.species)
        exit(0)

    if args.info:
        summary = 'Species name %s:\n' % args.species 
        for factor, _list in genomes[args.species].items():
            summary += '{}:\t'.format(factor)
            for i in range(len(indecies)):
                summary += base[i] + ': ' + _list[indecies[i]] + '\t'
            summary += '\n'
        print(summary)
        exit(0)

    for site in websites:
        cmd = [['wget']]
        prefix = site.split('/')[-1]

        if args.get_ann:
            nfile = prefix + '_genomic.gtf.gz'
        else:
            nfile = prefix + '_genomic.fna.gz'
            nfile_ix = prefix + '_assembly_report.txt'
            cmd.append(['wget', site + '/' + nfile_ix, '-P', args.out_dir])

        cmd[0] += [site + '/' + nfile, '-P', args.out_dir]
        try:
            for c in cmd:
                subprocess.call(c)
        except:
            print("Download Failed")
        else:
            print('Completed:{}'.format(nfile))
        
        if args.get_ann or args.dont_convert:
            exit(0)

        geno_fn = args.out_dir + nfile
        geno_ix_fn = args.out_dir + nfile_ix
        subprocess.call(['gunzip', geno_fn])

        convert_names(geno_fn.replace('.gz', ''), geno_ix_fn, args.use_ucsc)
