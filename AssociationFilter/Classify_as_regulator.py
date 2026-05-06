import os, sys, shutil, pickle, gzip, argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, fcluster
from itertools import combinations
import tempfile


BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Setup paths
target_dir = os.path.abspath('/home/pgohl/ModCRE_Package/exe')
sys.path.append(target_dir)
import fimo

# Constants
GENOME_DIR = "/home/pgohl/Work/pgohl/mutations/AdAstra/hg19/"
PWM_DIR = os.path.expanduser("~/FinalModCRElibTest/ModCRElib/ExternalPWMs/jaspar_pwms")

# --- Helper Functions ---

def load_tf_mapping(path=os.path.join(BASE_DIR, "TF_to_JASPAR.pkl")):
    with open(path, "rb") as fd:
        return pickle.load(fd)

def jaccard_similarity(list1, list2):
    set1, set2 = set(list1), set(list2)
    union = len(set1.union(set2))
    return len(set1.intersection(set2)) / union if union > 0 else 1.0


def load_chromosome_sequence(filepath):
    seq = []
    with gzip.open(filepath, "rt") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return ''.join(seq)


def extract_snp_windows(
    df,
    genome_dir,
    window=1000
):
    """
    df must contain columns: SNP, SNPChr, SNPPos
    Returns dict: {snp_id: (chrom, pos, sequence)}
    """
    snps_by_chr = defaultdict(list)
    for _, row in df.iterrows():
        chrom = f"chr{row['SNPChr']}"
        pos = int(row["SNPPos"])
        snps_by_chr[f"{chrom}.fa.gz"].append((row["SNP"], pos))
    results = {}
    for chrom_file, snplist in snps_by_chr.items():
        chrom_path = os.path.join(genome_dir, chrom_file)
        if not os.path.exists(chrom_path):
            print(f"WARNING: Missing genome file: {chrom_path}")
            continue
        seq = load_chromosome_sequence(chrom_path)
        chrom_len = len(seq)
        for snpid, pos in snplist:
            start = max(1, pos - window)
            end   = min(chrom_len, pos + window)
            subseq = seq[start - 1 : end]
            results[snpid] = (chrom_file[3:-6], pos, subseq.upper())
    return results

def run_fimo(snp_windows, window, TFs, TF_TO_JASPAR, workdir):
    known_regions_motifs = []

    for snpid, (chrom, pos, seq) in snp_windows.items():
        region_binds = []

        out_fasta = os.path.join(workdir, "input.fa")
        with open(out_fasta, "w") as out:
            out.write(f">{snpid}|chr{chrom}|pos{pos}|window{window}\n{seq}\n")

        pwm_tmp = os.path.join(workdir, "PWMs")
        os.makedirs(pwm_tmp, exist_ok=True)
        #pwm_tmp="."
        meme_db = os.path.join(pwm_tmp, "database.meme")
        with open(meme_db, "wb") as dest:
            for key in TFs:
                if key in TF_TO_JASPAR:
                    with open(f"{PWM_DIR}/{TF_TO_JASPAR[key]}.meme", "rb") as src:
                        shutil.copyfileobj(src, dest)
        #print(workdir,meme_db)
        #exit(0)
        x = fimo.get_fimo_obj(
            meme_db,
            out_fasta,
            fimo_pvalue_threshold=0.00001,
            dummy_dir=workdir
        )

        for hit in x.get_hits():
            region_binds.append(hit.get_hit())

        known_regions_motifs.append(region_binds)

    return known_regions_motifs


# --- Main Logic Wrapper ---

def analyze_regulatory_regions(gene, chrom, tf_list, query_tf, query_pos: int, eqtl_file):
    """
    Core logic to be called by other scripts.
    Returns dictionary of clustering and similarity results.
    """
    TF_TO_JASPAR = load_tf_mapping()
    with tempfile.TemporaryDirectory(prefix="assoc_") as workdir:

        # 1. Load eQTLs
        df = pd.read_csv(eqtl_file, sep='\t', usecols=["SNP", "SNPChr", "SNPPos", "GeneSymbol"])
        df = df[(df["GeneSymbol"] == gene) & (df["SNPChr"] == chrom)]
    
        if df.empty:
            return {"error": f"No eQTLs found for {gene} on chrom {chrom}"}

        # 2. Extract Motifs for all eQTL regions
        results_map = extract_snp_windows(df, GENOME_DIR)
        # No filtering: all motifs returned are kept
        known_motifs = run_fimo(results_map, 1000, tf_list, TF_TO_JASPAR, workdir)

        # 3. Clustering (Ward linkage + Max Clusters to prevent sparse results)
        sets = [set(m) for m in known_motifs]
        n = len(sets)
        if n < 2:
            return {"error": "Not enough regions to perform clustering."}

        dist_mat = np.zeros((n, n))
        for i, j in combinations(range(n), 2):
            dist_mat[i, j] = dist_mat[j, i] = 1 - jaccard_similarity(sets[i], sets[j])
    
        # Ward minimizes within-cluster variance
        linked = linkage(dist_mat[np.triu_indices(n, 1)], method='ward')
    
        # Force 5 clusters to ensure population density; adjust as needed
        cluster_ids = fcluster(linked, 5, criterion='maxclust')

        groups = defaultdict(list)
        for i, cid in enumerate(cluster_ids):
            groups[cid].append(known_motifs[i])

        # 4. Process the New (Query) Region
        asb_df = pd.DataFrame([["ASB", chrom, query_pos, gene]], columns=['SNP','SNPChr','SNPPos','GeneSymbol'])
        new_res = extract_snp_windows(asb_df, GENOME_DIR)
        new_motifs_list = run_fimo(new_res, 1000, tf_list, TF_TO_JASPAR, workdir)
    
        if not new_motifs_list:
            return {"error": "No motifs found in the query region."}
    
        new_set = set(new_motifs_list[0])

        # 5. Classification (Compare new region to cluster consensus)
        max_sim = -1
        best_c = None
    
        for cid, m_lists in groups.items():
            consensus = set().union(*m_lists)
            sim = jaccard_similarity(new_set, consensus)
            if sim > max_sim:
                max_sim, best_c = sim, cid

        return {
            "best_cluster_id": int(best_c),
            "similarity_score": round(float(max_sim), 4),
            "n_regions_in_cluster": len(groups[best_c]),
            "total_eqtl_regions": n,
            "is_regulator": max_sim >= 0.5,
            "cluster_map": {int(k): len(v) for k, v in groups.items()}
        }

# --- Standard execution for Argparse ---
if __name__ == "__main__":
    # [Argparse setup remains exactly as you have it]
    parser = argparse.ArgumentParser()
    # ... args ...
    args = parser.parse_args()
    
    tfs = args.tf.split(",")
    output = analyze_regulatory_regions(
        args.gene, args.chrom, tfs, args.queryTF, args.position, 
        os.path.join(BASE_DIR, '2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt')
    )
    print(output)
