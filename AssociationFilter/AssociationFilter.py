import requests
import pandas as pd
import pickle
#from Classify_as_regulator import analyze_regulatory_regions
import os
from pyliftover import LiftOver

BASE_DIR = os.path.dirname(os.path.abspath(__file__))



# Load the file
ABCdf = pd.read_csv(
    "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt",
    sep="\t"
)

with open(os.path.join(BASE_DIR, "TF_to_JASPAR.pkl"),"rb") as fd:
    TF_TO_JASPAR=pickle.load(fd)


def map_targets_by_tf(file_path, target_tf):
    """
    Searches a regulatory file for a specific TF (by UniprotID or NCBI GeneID)
    and maps all associated UniprotID.Target values to a dictionary.

    Args:
        file_path (str): Path to the regulatory data file (TSV/tab-separated).
        target_tf (str): The Transcription Factor ID (e.g., 'Q9H9S0' or '79923') to search for.

    Returns:
        dict: A dictionary where the key is the target_tf and the value is 
              a list of unique UniprotID.Target values.
    """
    # 1. Read the data (assuming the data is tab-separated)
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return {}        
    # Standardize column names by stripping extra whitespace if present
    df.columns = df.columns.str.strip()
    # Define the columns to search
    tf_search_cols = ['UniprotID.TF', 'Name.TF']
    # Ensure NCBI Gene ID is treated as string for exact match
    df['Name.TF'] = df['Name.TF'].astype(str)
    # 2. Create the search condition (Boolean Mask)
    # Check if the target_tf is in either UniprotID.TF OR NCBI.GeneID.TF column
    condition = df[tf_search_cols].eq(target_tf).any(axis=1)
    # 3. Filter the DataFrame
    filtered_df = df[condition]
    # 4. Extract and store the results in the dictionary format
    if not filtered_df.empty:
        # Get the unique list of UniprotID.Target values
        #target_list = filtered_df['UniprotID.Target'].astype(str).unique().tolist()
        target_list = filtered_df['Name.Target'].astype(str).unique().tolist()
        # Store in the required dictionary format
        results = {target_tf: target_list}
    else:
        results = {target_tf: []}
    return results


def map_tfs_by_target(file_path, target_gene):
    """
    Searches a regulatory file for a specific target gene (by GeneSymbol or Ensembl ID)
    and maps all associated TFs to a dictionary.

    Args:
        file_path (str): Path to the regulatory data file (TSV/tab-separated).
        target_gene (str): Target gene identifier (e.g. 'CLEC12A' or 'ENSG00000172322').

    Returns:
        dict: A dictionary where the key is the target_gene and the value is
              a list of unique TF names.
    """
    # 1. Read the data
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return {}

    # Standardize column names
    df.columns = df.columns.str.strip()
    # Columns to search for target gene
    target_search_cols = ['UniprotID.Target', 'Name.Target']
    # Ensure string comparison
    df['Name.Target'] = df['Name.Target'].astype(str)
    # 2. Create Boolean mask
    condition = df[target_search_cols].eq(target_gene).any(axis=1)
    # 3. Filter DataFrame
    filtered_df = df[condition]
    # 4. Extract TFs
    if not filtered_df.empty:
        tf_list = filtered_df['Name.TF'].astype(str).unique().tolist()
        results = {target_gene: tf_list}
    else:
        results = {target_gene: []}
    return results


def get_genes_near_snp(chrom, pos, window=500000, assembly='hg38'):
    """
    Query Ensembl REST API for genes whose TSS is within ±window of SNP.
    Returns list of (gene_id, gene_name, tss, strand, distance).
    """
    if assembly == 'hg19':
        lo = LiftOver('hg19', 'hg38')
        coord = lo.convert_coordinate(chrom, int(pos))
        pos = coord[0][1]
    
    server = "https://rest.ensembl.org"
    region_start = pos - window
    region_end = pos + window
    ext = f"/overlap/region/human/{chrom}:{region_start}-{region_end}?feature=gene"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    genes = r.json()
    hits = []
    for g in genes:
        gene_id = g["id"]  # Ensembl stable gene ID, e.g. ENSG00000123456
        gene_name = g.get("external_name", gene_id)  # Symbol or fallback to ID
        strand = g["strand"]
        start = g["start"]
        end = g["end"]
        # Determine TSS
        tss = start if strand == 1 else end
        # Check TSS within ±window
        if abs(tss - pos) <= window:
            hits.append((gene_id, gene_name, tss, strand, abs(tss - pos)))
    # Sort by distance to SNP
    hits.sort(key=lambda x: x[4])
    return hits

def gene_distance_score(distance_bp):
    if pd.isna(distance_bp) or distance_bp > 500_000:
        return 0.0
    bins = [0, 100_000, 200_000, 300_000, 400_000, 500_000]
    scores = [1.0, 0.8, 0.6, 0.4, 0.2]
    for i in range(len(bins) - 1):
        if bins[i] <= distance_bp < bins[i + 1]:
            return scores[i]
    return 0.0

def find_target_genes(chrom, coord):
    matches = ABCdf[
        (ABCdf["chr"] == chrom) &
        (ABCdf["start"] < coord) &
        (ABCdf["end"] > coord)
    ]

    return matches["TargetGene"].unique()



def run_regulatory_analysis(
    TF,
    chrom,
    SNPeBoT2_prediction,
    pos,
    assembly='hg38',
    distancerange=10000,
    file_path_TFlink=os.path.join(BASE_DIR,"TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv"),
    disgenet_file=os.path.join(BASE_DIR,"disgenet/all_gene_disease_pmid_associations.tsv.gz")
):
    """
    Main analysis entry point.
    Returns a pandas DataFrame (final_df).
    """

    start = pos - distancerange
    end = pos + distancerange

    # --- Load DisGeNET once ---
    disgenet_df = pd.read_csv(
        disgenet_file,
        sep="\t",
        compression="gzip",
        usecols=["geneSymbol", "diseaseName", "diseaseId", "pmid"]
    )
    if chrom != "X" and chrom != "Y" and chrom != "M":
        chrom = int(chrom.replace("chr",""))
    # --- TF → target genes ---
    target_mapping_name = map_targets_by_tf(file_path_TFlink, TF)
    target_genes = set(target_mapping_name.get(TF, []))

    if not target_genes:
        print("no target genes for the given TF")
        return pd.DataFrame()  # early exit

    # --- Nearby genes ---
    genes = get_genes_near_snp(chrom, pos, 500000, assembly)
    if assembly == 'hg38':
        lo = LiftOver('hg38', 'hg19')
        coord = lo.convert_coordinate('chr'+str(chrom), int(pos))
        hg19pos = coord[0][1]
    else:
        hg19pos = pos

    regulated_genes = find_target_genes('chr'+str(chrom), hg19pos)


    rows = []

    for gene_id, gene_name, tss, strand, dist in genes:

        if gene_name not in target_genes:
            continue

        within_distance = dist <= distancerange

        # --- CREM / regulatory logic ---
        TFList = map_tfs_by_target(file_path_TFlink, gene_name)
        if gene_name in regulated_genes:
            is_regulator = True
            regulator_score = 1
        else:
            is_regulator = False
            regulator_score = 0
        

        disease_hits = disgenet_df[disgenet_df["geneSymbol"] == gene_name]

        if disease_hits.empty:
            rows.append({
                "TF": TF,
                "TargetGeneSymbol": gene_name,
                "TargetGeneID": gene_id,
                "SNPeBoT2_prediction": SNPeBoT2_prediction,
                "gene_distance_bp": dist,
                "regulator_score": regulator_score,
                "regulator": is_regulator,
                "disease_name": None,
                "disease_id": None,
                "pmid": None
            })
        else:
            for _, d in disease_hits.iterrows():
                rows.append({
                    "TF": TF,
                    "TargetGeneSymbol": gene_name,
                    "TargetGeneID": gene_id,
                    "SNPeBoT2_prediction": SNPeBoT2_prediction,
                    "gene_distance_bp": dist,
                    "regulator_score": regulator_score,
                    "regulator": is_regulator,
                    "disease_name": d["diseaseName"],
                    "disease_id": d["diseaseId"],
                    "pmid": d["pmid"]
                })

    final_df = pd.DataFrame(rows)

    if final_df.empty:
        print("there is an empty final df")
        return final_df

    # --- Merge PMIDs ---
    group_cols = [
        "TargetGeneSymbol",
        "TargetGeneID",
        "SNPeBoT2_prediction",
        "gene_distance_bp",
        "regulator_score",
        "regulator",
        "disease_name",
        "disease_id",
    ]

    final_df = (
        final_df
        .groupby(group_cols, dropna=False, as_index=False)
        .agg({
            "pmid": lambda x: ";".join(sorted(set(x.dropna().astype(str))))
        })
    )

    # --- Completeness score ---
    final_df["completeness_score"] = (
        final_df["SNPeBoT2_prediction"].isin(["gain", "loss"]).astype(int)
        + final_df["regulator_score"]
        + final_df["regulator"].astype(int)
        + final_df["disease_name"].notna().astype(int)
    )

    final_df = final_df.sort_values(
        by=[
            "completeness_score",
            #"gene_distance_bp",
            "regulator",
            "disease_name"
        ],
        ascending=[False, False, True]
    ).reset_index(drop=True)

    return final_df



import sys

if __name__ == "__main__":
    TF = sys.argv[1]
    chrom = sys.argv[2]
    prediction = sys.argv[3]
    pos = int(sys.argv[4])
    assembly = sys.argv[5]
    output_file = sys.argv[6]
    print(output_file)

    df = run_regulatory_analysis(
        TF=TF,
        chrom=chrom,
        SNPeBoT2_prediction=prediction,
        pos=pos,
        distancerange=500000
    )

    df.to_csv(output_file, index=False)


