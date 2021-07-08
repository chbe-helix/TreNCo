import pandas as pd
import glob
import numpy as np
import os
import multiprocessing


mean_gene_exp = None
mean_enh_k27ac = None
context_info = ''



def init(_exp, _enh, _context):
    global mean_gene_exp
    global mean_enh_k27ac
    global context_info
    mean_gene_exp = _exp
    mean_enh_k27ac = _enh
    context_info = _context

def get_sample_exp_df(exp_df, tissue, time, excl):
    exp_ft_df = exp_df.ix[:,[all(sb not in c for sb in excl) and (tissue in c) and (time in c) for c in exp_df.columns]]
    mean_exp_ft_df = pd.DataFrame(exp_ft_df.mean(axis=1))
    #print("Expression:", exp_ft_df.shape)
    #print(exp_ft_df.columns)
    return mean_exp_ft_df

def get_sample_k27ac_df(k27ac_df, tissue, time):
    k27ac_ft_df = k27ac_df.ix[:,[(tissue in c) and (time in c) for c in k27ac_df.columns]]
    mean_k27ac_ft_df = pd.DataFrame(k27ac_ft_df.mean(axis=1))
    #print("K27ac:", k27ac_ft_df.shape)
    #print(k27ac_ft_df.columns)
    return mean_k27ac_ft_df

def get_distance_df(gene_df,enh_df):
    enh_pos = []
    for r in enh_df.iterrows():
        pos = "{}:{}-{}".format(r[1][0],r[1][1],r[1][2])
        enh_pos.append(pos)
    gene_names = list(gene_df[3])
    out_d = pd.DataFrame(0,index=enh_pos,columns=gene_names)
    for g in gene_df.iterrows():
        name = g[1][3]
        if g[1][5] == '-':
            g_start = g[1][2]
        else:
            g_start = g[1][1]
        for e in enh_df.iterrows():
            e_pos = "{}:{}-{}".format(e[1][0],e[1][1],e[1][2])
            e_mid = e[1][1]+((e[1][2]-e[1][1])/2)
            d = np.abs(g_start - e_mid)
            out_d.ix[e_pos,name] = int(d)
    return out_d

def get_distance_based_weight(distance):
    ##set weight 0 for distance greater than 500 kb
    distance[distance > 500000] = 1e10000
    dist_weight = 1/np.log2(distance)
    if set(dist_weight.values.flatten()) == {0}:
        wmin = np.sort(list(set(dist_weight.values.flatten())))[0]
        wmax = np.sort(list(set(dist_weight.values.flatten())))[0]
        scaled_dist_weight = dist_weight
    else:
        wmin = np.sort(list(set(dist_weight.values.flatten())))[1]
        wmax = np.sort(list(set(dist_weight.values.flatten())))[-1]
        scaled_dist_weight = (dist_weight-wmin)/(wmax-wmin)
        scaled_dist_weight[scaled_dist_weight<0] = 0
    return(scaled_dist_weight)

def get_enh_gene_weights(tup):
    tad_dir, context_info = tup
    enh_path = '{}/enh.bed'.format(tad_dir)
    gene_path = '{}/gene.bed'.format(tad_dir)

    if os.path.isfile(enh_path) and os.path.isfile(gene_path):
        enh_df = pd.read_csv(enh_path,sep='\t',header=None)
        gene_df = pd.read_csv(gene_path,sep='\t',header=None)
        enh_gene_dist = get_distance_df(gene_df,enh_df)
        scaled_weight = get_distance_based_weight(enh_gene_dist)
        #scaled_weight.to_csv('{}/enh_gene_dist_scaled_weight.txt'.format(tad_dir),sep='\t')
#         m_gene_exp = share_data.g_exp
#         m_enh_k27ac = share_data.e_sig
        tad_gene_exp = mean_gene_exp.ix[scaled_weight.columns,:]
        tad_enh_signal = mean_enh_k27ac.ix[scaled_weight.index,:]
        enh_gene_weight = np.multiply(
            scaled_weight,np.sqrt(np.matmul(tad_enh_signal,tad_gene_exp.T)))
        enh_gene_weight.to_csv('{}/{}_enh_gene_weight.txt'.format(tad_dir,context_info),sep='\t')

def fill_network(d):
    enh_path = '{}/enh.bed'.format(d)
    gene_path = '{}/gene.bed'.format(d)
    weight_path = '{}/{}_enh_gene_weight_minus_expected_tanh.txt'.format(d, context_info)
    if os.path.isfile(enh_path) and os.path.isfile(gene_path) and os.path.isfile(weight_path):
        weight = pd.read_csv(weight_path, sep='\t',index_col=0)
        # weight_dict = weight.to_dict()
        # for c in weight_dict.keys():
        #     for r in weight_dict[c].keys():
        #         all_network.ix[r,c] = np.float32(weight_dict[c][r])
        return weight
    else:
        return None

if __name__ == "__main__":

    data_analysis = {
        'embryonic-facial-prominence':
            ['e11.5', 'e15.5', 'e12.5', 'e13.5', 'e14.5'],
        'forebrain':
            ['e12.5', 'e11.5', 'e15.5', 'e13.5', 'e14.5', 'e16.5', 'P0'],
        'heart':
            ['e11.5', 'e16.5', 'e13.5', 'e15.5', 'P0', 'e12.5', 'e14.5'],
        'hindbrain':
            ['e14.5', 'e11.5', 'e15.5', 'e13.5', 'P0', 'e12.5', 'e16.5'],
        'intestine':
            ['P0', 'e15.5', 'e14.5', 'e16.5'],
        'kidney':
            ['e16.5', 'e14.5', 'e15.5', 'P0'],
        'limb':
            ['e12.5', 'e13.5', 'e15.5', 'e14.5', 'e11.5'],
        'liver':
            ['P0', 'e16.5', 'e11.5', 'e15.5', 'e12.5', 'e13.5', 'e14.5'],
        'lung':
            ['e14.5', 'P0', 'e16.5', 'e15.5'],
        'midbrain':
            ['e16.5', 'e13.5', 'e14.5', 'e11.5', 'e15.5', 'e12.5', 'P0'],
        'neural-tube':
            ['e14.5', 'e13.5', 'e11.5', 'e15.5', 'e12.5'],
        'stomach':
            ['e14.5', 'P0', 'e15.5', 'e16.5']}

    exp_pc_df = pd.read_csv('/home/shared/Data/encode/mouse/gene_expression_log2TPM_signal_processed_matrix.txt',sep='\t',index_col=0)
    enh = pd.read_csv('/home/shared/Data/encode/mouse/enh_H3K27ac_log2TPM_signal_processed_matrix.txt',sep='\t',index_col=0)

    ## pariwise pearson correlation coefficient amoung its techincal and biological replicate was < 0.3. Therefore going to remove these samples
    rnaseq_to_remove = ['ENCFF662WLV','ENCFF705YYN','ENCFF026IQF','ENCFF892XES','ENCFF668RJN','ENCFF993NNK','ENCFF569ODO','ENCFF863BCB','ENCFF250DXJ','ENCFF434WEQ']
    k27ac_to_remove = ['ENCFF001KKQ','ENCFF001KKP','ENCFF001KME','ENCFF001KMD','ENCFF001KGS','ENCFF001KGT']


    tad_dirs = glob.glob('/home/shared/Data/encode/mouse/heart/TADs/TAD_*')
    count = 0
    #take one network index and sort all other networks using same index
    hdf = pd.HDFStore('/home/shared/Data/encode/mouse/heart/old/heart_P0_enh-gene_network.h5', mode='r')
    all_weight = hdf['/enh_gene_all_weights']
    hdf.close()
    out_index = all_weight.index
    out_columns = all_weight.columns

    del(all_weight)
    # for tissue in data_analysis.keys():
    #     for age in data_analysis[tissue]:
    #         context_info = '{}_{}'.format(tissue,age)
    #         if context_info in ['heart_e11.5','heart_e16.5','heart_e13.5']:
    #             next
    #         count += 1
    #         print(count, context_info)
    #         # take mean of all techincal and biological replicates
    #         mean_gene_exp = get_sample_exp_df(
    #             exp_pc_df, age, tissue, rnaseq_to_remove)
    #         mean_enh_k27ac = get_sample_k27ac_df(enh, age, tissue)
    #
    #         print('Calculating weights for each TAD')
    #         all_inputs = []
    #
    #         for d in tad_dirs:
    #             all_inputs.append((d,context_info))
    #
    #         p = multiprocessing.Pool(64, initializer = init, initargs = (mean_gene_exp, mean_enh_k27ac, context_info))
    #
    #         out = p.map(get_enh_gene_weights, all_inputs)
    #
    #         del(out)
    #
    #
    print('Calculate weights of Enh-Gene by subtracting expected')

    for file in tad_dirs:

        tad_files = glob.glob('{}/*_enh_gene_weight.txt'.format(file))
        for indx,c in enumerate(tad_files):
            if indx == 0:
                sum_enh_gene = pd.read_csv(c,sep='\t',index_col=0)

            else:
                other_enh_gene = pd.read_csv(c, sep='\t',index_col=0)
                sum_enh_gene += other_enh_gene

        if len(tad_files) > 66:
            print('Inconsistent number of context files. The number of files are greater than 66')
            sys.exit()

        expected = sum_enh_gene/len(tad_files)

        for c in tad_files:
            c_enh_gene = pd.read_csv(c,sep='\t',index_col=0)
            w_enh_gene = np.tanh(c_enh_gene - expected)
            path = c.split('.txt')[0]
            w_enh_gene.to_csv('{}_minus_expected_tanh.txt'.format(path),sep='\t')


    for tissue in data_analysis.keys():
        for age in data_analysis[tissue]:
            context_info = '{}_{}'.format(tissue, age)
            count += 1
            print(count, context_info)
            #init(None, None, context_info)
            p = multiprocessing.Pool(64, initializer = init, initargs = (None, None, context_info))

            dfs = p.map(fill_network, tad_dirs)

            print('All DFS', len(dfs))
            p.close()
            print('Gathering all weights from all DFs')
            all_network = pd.concat(dfs)

            all_network.fillna(0, inplace=True)
            all_network = all_network.ix[out_index, out_columns]
            # #print('Sum of all weights', np.sum(np.array(all_network)))
            hdf = pd.HDFStore('/home/shared/Data/encode/mouse/heart/{}_enh-gene_network_minus_expected_tanh_final.h5'.format(context_info))
            hdf['/enh_gene_all_weights'] = all_network
            hdf.close()

            print(all_network.shape)
            print('======================================')
