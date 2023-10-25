# sc
import pegasus as pg
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import scdrs

# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# data
import numpy as np
import pandas as pd
import re
from io import StringIO
import mygene

def get_ens_dict(file_path):
    with open(file_path) as f:
        gtf = list(f)
    gtf = [x for x in gtf if not x.startswith('#')]
    gtf = [x for x in gtf if 'gene_id "' in x and 'gene_name "' in x]
    if len(gtf) == 0:
        print('you need to change gene_id " and gene_name " formats')
    gtf = list(map(lambda x: (x.split('gene_id "')[1].split('"')[0], x.split('gene_name "')[1].split('"')[0]), gtf))
    gtf = dict(set(gtf))
    return gtf

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def replace(my_list, my_dict):
    return [x if x not in my_dict else my_dict[x] for x in my_list]

### load h5ad dataset
data=sc.read_h5ad("/path/to/data.h5ad")
cases=data[data.obs.dx=='AD']

### get gene symbol, name, entrez ID
gtf_dict = get_ens_dict('/path/to/Homo_sapiens.GRCh38.107.chr.gtf') #replace with your file path
genes=list(gtf_dict.keys())
conversion_table=mygene.MyGeneInfo().getgenes(genes, fields='name,symbol,entrezgene,taxid', as_dataframe=True)
dict_genes=pd.Series(conversion_table['symbol'].values,index=conversion_table['_id']).to_dict()

### Creating standard geneset as scDRS input using MAGMA outputs
numbers = re.compile(r'(\d+)')
files = ['ensemblProtCodGenes35kb10kbNoChryNoBmhc__edu.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__neu.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__height.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__drinking.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__bmi.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__dm2.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__cad.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__uc.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__cd.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__ibd.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__sle.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__ra.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__ds.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__sz3.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__mdd.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__bip.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__asd.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__ms.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__pd.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__als.genes.out','ensemblProtCodGenes35kb10kbNoChryNoBmhc__alzBellenguezNoApoe.genes.out']
panda_df=[]
files_df=[]

for f in files:
    df = pd.read_csv(f, sep='\s+')
    df=df.sort_values('ZSTAT', ascending=False). head(1000)
    df_cut=df[['GENE', 'ZSTAT']]
    repl =replace(df_cut.GENE,gtf_dict)
    df_cut['GENE'] = repl
    df_cut=df_cut[~df_cut.GENE.str.startswith('ENSG')]
    mylist=[]
    for key, value in zip(df_cut.GENE, df_cut.ZSTAT):
        mylist.append(key)
        mylist.append(':')
        mylist.append(value)
        mylist.append(',')
    s = StringIO()
    print(*mylist, sep="", file=s)
    panda_df.append(s.getvalue())
    f=f.replace('ensemblProtCodGenes35kb10kbNoChryNoBmhc__','')
    f=f.replace('BellenguezNoApoe.genes.out','')
    files_df.append(f.replace('.genes.out',''))
                    
genes=pd.DataFrame(panda_df)
exper=pd.DataFrame(files_df)
final=pd.concat([exper, genes], axis=1)
final.columns=['exp','genes']
final['genes']=final['genes'].str.replace(r'.$', '')
final.to_csv("standard_geneset.gs", sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE, escapechar=' ')

### Performing scDRS analysis on cases
dict_gs = scdrs.util.load_gs('/path/to/standard_geneset.gs', src_species="human", dst_species="human", to_intersect=cases.var_names)
dict_gs={k:v for k,v in dict_gs.items()}

scdrs.preprocess(cases, n_mean_bin=20, n_var_bin=20, copy=False)

dict_df_score = dict()
for trait in dict_gs:
    gene_list, gene_weights = dict_gs[trait]
    dict_df_score[trait] = scdrs.score_cell(data=cases, gene_list=gene_list, gene_weight=gene_weights, ctrl_match_key="mean_var", n_ctrl=100, weight_opt="vs", return_ctrl_raw_score=False, return_ctrl_norm_score=True, verbose=False)
    dict_df_score[trait].to_csv("/path/to/MAGMA/{}_microglia_cases_subtype_standardgeneset.csv".format(trait), sep="\t", quoting=csv.QUOTE_NONE,  escapechar=' ')
	df_stats = scdrs.method.downstream_group_analysis(adata=cases, df_full_score=dict_df_score[trait], group_cols=["subtype"])["subtype"]
    display(df_stats.style.set_caption("Group-level statistics for {}".format(trait)))
    df_stats.to_csv("/path/to/MAGMA/{}_microglia_cases_df_stats_subtype_standardgeneset.csv".format(trait), sep="\t", quoting=csv.QUOTE_NONE,  escapechar=' ')
    dict_microg=pd.Series(dict_df_score[trait]['zscore'],index=dict_df_score[trait].index).to_dict()
    cases.obs['scdrs_zscore'] = cases.obs.index.map(dict_microg)
    cases.obsm['umap']=cases.obsm['X_umap']
    sc.pl.umap(cases, color='scdrs_zscore',cmap = sns.color_palette("cubehelix", as_cmap=True))