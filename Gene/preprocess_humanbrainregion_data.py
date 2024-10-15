#!/usr/bin/env python

"""
convert h5 format (downloaded) to h5ad
(in babel conda environment)
"""

import h5py
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix
import random
import numpy as np

#############################################################################
## load and preprocess data
#############################################################################

data_dir = '/net/noble/vol8/ranz0/aneuploidy/data/'
## read the data
data = h5py.File(data_dir + 'linnarsson-lab/HumanFetalBrainPool.h5', 'r')

## output gene and cell annotations
df_cell = pd.DataFrame({'Doner':data['shoji']['Donor'][:], 'Age':data['shoji']['Age'][:], 'Sex':data['shoji']['Sex'][:], 'batch':data['shoji']['Chemistry'][:], 'Region':data['shoji']['Region'][:], 'Subregion':data['shoji']['Subregion'][:], 'Tissue':data['shoji']['Tissue'][:], 'CellClass':data['shoji']['CellClass'][:], 'CellID':data['shoji']['CellID'][:]})
df_cell.index = df_cell.CellID
df_cell.index.name = 'index'
df_cell.to_csv(data_dir + 'linnarsson-lab/cell_annot.txt', sep='\t')

df_gene = pd.DataFrame({'Gene':data['shoji']['Gene'][:], 'Accession':data['shoji']['Accession'][:], 'Chromosome':data['shoji']['Chromosome'][:], 'Start':data['shoji']['Start'][:], 'End':data['shoji']['End'][:], 'ValidGenes':data['shoji']['ValidGenes'][:], 'SelectedFeatures':data['shoji']['SelectedFeatures'][:]})
df_gene.index = df_gene.Accession
df_gene.index.name = 'index'
df_gene.to_csv(data_dir + 'linnarsson-lab/gene_annot.txt', sep='\t')

## filter by genes in our dataset:
rna_data_gene = pd.read_csv(data_dir +'Disteche_RNA3-042-nova_data/matrices/115.020b_gene_annotations.txt', 
                            delimiter='\t', header=None)
rna_data_gene.columns = ['Accession','Gene']


#############################################################################
## save to h5ad format to be input to scanpy
#############################################################################
"""
it turns out the matrix is huge and I don't have a good way to directly convert it, so instead I returned two versions: 
1. subset of cells for preliminary exploration; 
2. store each chunk of data separately
"""

## 1. subsample nsubsample cells from the whole matrix, for preliminary exploration
nsubsample = 200000
random.seed(101)
selected_index = sorted(random.sample(list(range(data['shoji']['Expression'].shape[0])), nsubsample))
testt = data['shoji']['Expression'][np.array(selected_index), :]
adata = ad.AnnData(csr_matrix(testt))
adata.obs = df_cell.iloc[selected_index]
adata.var = df_gene
adata = adata[:, adata.var['Gene'].isin(rna_data_gene['Gene'])]
adata.write(filename = data_dir + 'linnarsson-lab/HumanFetalBrainPool_sub' + str(nsubsample)+'.h5ad', compression=None, compression_opts=None, force_dense=None)


## 2. save the whole Expression matrix as h5ad
# the problem is h5 file is so big, so we need to convert them per nk cells, save and concatenate later
adatas = []
nk = 200000
nfile = data['shoji']['Expression'].shape[0] // nk + 1
for i in range(nfile):
	print(i)
	testt = data['shoji']['Expression'][i*nk : min((i+1)*nk, data['shoji']['Expression'].shape[0]), :]
	adata = ad.AnnData(csr_matrix(testt))
	adata.obs = df_cell.iloc[i*nk : min((i+1)*nk, data['shoji']['Expression'].shape[0])]
	adata.var = df_gene
	adata.write(filename = data_dir + 'linnarsson-lab/HumanFetalBrainPool' + str(i)+'.h5ad', compression=None, compression_opts=None, force_dense=None)
	#adatas.append(adata)
	del testt
	del adata

del data
del df_cell
del df_gene

# concatenate all of them - TODO
adatas = []
nfile = 9
for i in range(nfile):
	adatas.append(ad.read_h5ad(data_dir + 'linnarsson-lab/HumanFetalBrainPool' + str(i)+'.h5ad'))

adata = adatas[0].concatenate(adatas[1::], )
adata.write(filename = data_dir + 'linnarsson-lab/HumanFetalBrainPool.h5ad', compression=None, compression_opts=None, force_dense=None)
