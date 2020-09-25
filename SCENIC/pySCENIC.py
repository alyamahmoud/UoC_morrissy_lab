import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar
from dask.distributed import Client
client = Client(processes=False)

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns
import json

DATA_FOLDER="/work/morrissylab/alyaa.mohamed/for_dmitri"
RESOURCES_FOLDER="/work/morrissylab/alyaa.mohamed/for_dmitri"
DATABASE_FOLDER = "/work/morrissylab/alyaa.mohamed/for_dmitri"
# SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.mc9nr.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_mgi_tfs.txt')
SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "exprMat.txt")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")

ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
ex_matrix.shape

MOTIFS_MGI_FNAME = os.path.join(DATA_FOLDER, 'motifs-v9-nr.mgi-m0.001-o0.0.tbl')
OUT_TFS_MGI_FNAME = os.path.join(DATA_FOLDER, 'mm_mgi_tfs.txt')
df_motifs_mgi = pd.read_csv(MOTIFS_MGI_FNAME, sep='\t')
mm_tfs = df_motifs_mgi.gene_name.unique()
with open(OUT_TFS_MGI_FNAME, 'wt') as f:
    f.write('\n'.join(mm_tfs) + '\n')
    
len(mm_tfs)

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs

tf_names = load_tf_names(MM_TFS_FNAME)

adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

with open ('adjacencies_test.json', 'w') as f:
	json.dump (adjacencies, f)

