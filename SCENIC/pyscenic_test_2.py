#! /usr/bin/env python
import os, sys, time

# =============================================================================================
if __name__ == '__main__':
        t0 = time.time()
        # =============================================================================================
        from dask.distributed import LocalCluster, Client
        print("Creating a cluster...")
        nproc = 40
        c = LocalCluster(processes=True, n_workers=0)
        print("Done.")
        print(c)
        print("Adding more workers...")
        for i in range(nproc):
                print("%d " % (i+1), end="")
                c.start_worker(ncores=1, local_dir="/tmp/dask-worker")
        print("")
        print("Done.")
        print(c)
        print("Creating a client...")
        client = Client(c)
        print("Done.")
        print(client)

        # =============================================================================================
        t1 = time.time()
        print("Loading modules...")

        import glob
        import pandas as pd
        import numpy as np
        import pickle

        from arboreto.utils import load_tf_names
        from arboreto.algo import grnboost2

        from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
        from pyscenic.utils import modules_from_adjacencies, load_motifs
        from pyscenic.prune import prune2df, df2regulons
        from pyscenic.aucell import aucell
        import seaborn as sns

        print("All modules are loaded.")
        print("")

        # =============================================================================================
        t2 = time.time()

        data_dir = "data"
        res_dir = "resources"
        db_dir = "db_dir"

        dbs_mask = os.path.join(db_dir, "mm9-*.mc9nr.feather")

        motifs_annotation_fname = os.path.join(res_dir, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
        mm_fname = os.path.join(res_dir, 'mm_mgi_tfs.txt')
        expr_fname = os.path.join(data_dir, "exprMat.txt")
        motifs_fname = os.path.join(data_dir, 'motifs.csv')
        regulons_fname = os.path.join(data_dir, 'regulons.p')
        modules_fname = os.path.join(data_dir, 'modules.p')

        # =============================================================================================
        ex_matrix = pd.read_csv(expr_fname, sep='\t', header=0, index_col=0).T
        print("Expr Shape: ", ex_matrix.shape)

        df_motifs_mgi = pd.read_csv(motifs_annotation_fname, sep='\t')
        print("Motifs Shape:", df_motifs_mgi.shape)

        mm_tfs = df_motifs_mgi.gene_name.unique()
        print("Unique genes:", len(mm_tfs))

        print("Writing genes to '%s' file..." % mm_fname)
        open(mm_fname, 'wt').write('\n'.join(mm_tfs) + '\n')
        print("Done.")


        # =============================================================================================
        def name(fname):
            return os.path.splitext(os.path.basename(fname))[0]
        # =============================================================================================

        db_fnames = glob.glob(dbs_mask)
        print("DB files:", db_fnames)

        dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
        print(dbs)

        tf_names = load_tf_names(mm_fname)
        print("Loaded %d tf_names." % len(tf_names))

        print("Computing adjacencies...")
        adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True, client_or_address=client)
        print("Done.")

        print(adjacencies)

        # =============================================================================================
        print("Closing client...")
        client.close()
        print("Done.")

        print("Closing cluster...")
        c.close()
        print("Done.")
        # =============================================================================================
        df_fname = "adjacencies_test.dataframe"

        print("Writing adjacencies to DF '%s' file..." % df_fname)
        print('Writing adjacencies')
        try:
                adjacencies.to_pickle(df_fname)
                print("Done.")
        except:
                print("Failed.")
        #adjacencies.to_csv(adjacencies, index=False, sep='\t')
        with open (df_fname, 'wb') as f:
                pickle.dump(adjacencies, f)

        # =============================================================================================
        modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
        with open (modules_fname, 'wb') as f:
                pickle.dump(modules, f)

        df = prune2df(dbs, modules, motifs_annotation_fname)
        len(df)
        df.head()

        df.to_csv(os.path.join(res_dir, 'regulomes_test.csv'))
        with open (motifs_fname, 'wb') as f:
                pickle.dump(df, f)
        # py_df = pd.read_csv(os.path.join(RESOURCES_FOLDER, "regulomes_test.csv"), index_col=[0,1], header=[0,1], skipinitialspace=True)        

        py_df = pd.read_csv(os.path.join(res_dir, "regulomes_test.csv"), index_col=[0,1], header=[0,1], skipinitialspace=True)
        py_df[('Enrichment', 'TargetGenes')] = list(map(eval, py_df[('Enrichment', 'TargetGenes')]))

        regulons = df2regulons(py_df)
        len(regulons)
        
        
        with open ('regulons.p', 'wb') as f:
                pickle.dump(regulons, f)

        py_name2targets = {r.transcription_factor: set(r.genes) for r in regulons} # equivalent to 2.6_regulons_asIncidMat.txt
        regulons_asIncidMat = pd.DataFrame.from_dict(py_name2targets, orient = 'index')
        regulons_asIncidMat.to_csv(os.path.join(res_dir, 'regulons_asIncidMat.csv'))

        py_df.columns = py_df.columns.droplevel(0)
        py_df = py_df.reset_index()
        def getdb(s):
                for elem in eval(s):
                        if elem in dbs_mask:
                                return dbs_mask[elem]
                raise ValueError
        def getmodule(s):
                for elem in eval(s):
                        if elem not in dbs_mask:
                                return elem
                raise ValueError
        #py_df['DbId'] = py_df['Context'].apply(getdb)
        py_df['ModuleType'] = py_df['Context'].apply(getmodule)
        py_df.to_csv(os.path.join(res_dir, 'auc_scores.csv'))
        # =============================================================================================
        auc_mtx = aucell(ex_matrix, regulons, num_workers = 40)
        auc_mtx.to_csv(os.path.join(res_dir, 'auc_mtx.csv'))
        # =============================================================================================
        t10 = time.time()
        print("")
        print("All done (%.1f sec)." % (t10 - t0))
        print("")
        print("Cluster creation time: %.1f sec" % (t1 - t0))
        print("  Module loading time: %.1f sec" % (t2 - t1))
        print("     Computation time: %.1f sec" % (t10 - t2))
        print("")
        # =============================================================================================

