import os
import pandas as pd

# --- Configuration ---
# Define your input intersect files
high_fst_file = "intersect_high_fst_outliers.tsv"
private_sa_file = "intersect_privateSA_go.tsv"
private_aus_file = "intersect_privateAUS_go.tsv"

# --- NEW: Define paths to BOTH AF_Gap files ---
HIGH_FST_AF_GAPS_FILE = os.path.expanduser("~/msc/99_ortho/high_fst_outliers_af_gaps_with_freq.tsv")
ALL_PRIVATE_AF_GAPS_FILE = os.path.expanduser("~/msc/99_ortho/all_private_snps_af_gaps_with_freq.tsv") # The new file you just generated

# Output summary file name
OUTPUT_SUMMARY_FILE = "go_term_snp_overlap_confident_summary.tsv" # Overwriting previous version

# --- Step 0: Load and Combine ALL AF_Gap data ---
print(f"Loading AF_Gap data for High Fst SNPs from: {HIGH_FST_AF_GAPS_FILE}")
try:
    snp_af_data_high_fst = pd.read_csv(HIGH_FST_AF_GAPS_FILE, sep="\t")
    snp_af_data_high_fst['unique_snp_key'] = snp_af_data_high_fst['CHROM'].astype(str) + ':' + snp_af_data_high_fst['POS'].astype(str)
    print(f"Loaded {len(snp_af_data_high_fst)} High Fst SNPs with AF_Gap data.")
except FileNotFoundError:
    print(f"ERROR: High Fst AF_GAPS_FILE not found at {HIGH_FST_AF_GAPS_FILE}. Please ensure it exists.")
    exit()
except Exception as e:
    print(f"ERROR loading High Fst AF_GAPS_FILE: {e}")
    exit()

print(f"Loading AF_Gap data for ALL Private SNPs from: {ALL_PRIVATE_AF_GAPS_FILE}")
try:
    snp_af_data_private = pd.read_csv(ALL_PRIVATE_AF_GAPS_FILE, sep="\t")
    snp_af_data_private['unique_snp_key'] = snp_af_data_private['CHROM'].astype(str) + ':' + snp_af_data_private['POS'].astype(str)
    print(f"Loaded {len(snp_af_data_private)} All Private SNPs with AF_Gap data.")
except FileNotFoundError:
    print(f"ERROR: All Private AF_GAPS_FILE not found at {ALL_PRIVATE_AF_GAPS_FILE}. Please ensure it exists.")
    exit()
except Exception as e:
    print(f"ERROR loading All Private AF_GAPS_FILE: {e}")
    exit()

# Combine both AF_Gap dataframes.
# Use drop_duplicates to handle any SNPs that might be present in both lists (e.g., private AND high Fst),
# keeping the first occurrence (which would typically be from the high_fst_af_gaps_file if loaded first).
# In practice, AF_Gap value should be identical for such SNPs regardless of source.
snp_af_data_combined = pd.concat([snp_af_data_high_fst, snp_af_data_private], ignore_index=True).drop_duplicates(subset=['unique_snp_key'])
print(f"Combined AF_Gap data for {len(snp_af_data_combined)} unique SNPs.")


# --- Step 1: Load each intersect file and assign an initial label ---
dfs_by_label = {}

def load_and_label_file(filepath, label):
    if os.path.exists(filepath):
        df = pd.read_csv(filepath, sep="\t", header=None)
        df.columns = [
            "snp_chr", "snp_start", "snp_end", "snp_id",
            "gene_chr", "gene_start", "gene_end", "go_busco"
        ]
        df["initial_source"] = label
        df['unique_snp_key'] = df['snp_id'].apply(lambda x: ':'.join(x.split(':')[:2]))
        print(f"Loaded {len(df)} rows from {filepath} with label '{label}'")
        return df
    else:
        print(f"⚠️ Warning: {filepath} not found. Skipping.")
        return pd.DataFrame()

dfs_by_label["HighFST"] = load_and_label_file(high_fst_file, "HighFST")
dfs_by_label["PrivateSA"] = load_and_label_file(private_sa_file, "PrivateSA")
dfs_by_label["PrivateAUS"] = load_and_label_file(private_aus_file, "PrivateAUS")

# --- Step 2: Combine all DataFrames and merge with the COMPREHENSIVE AF_Gap data ---
all_combined_dfs = [df for df in dfs_by_label.values() if not df.empty]

if not all_combined_dfs:
    print("❌ No valid input files found for processing!")
    exit()

combined_raw = pd.concat(all_combined_dfs, ignore_index=True)

# Merge with the combined AF_Gap data
combined_with_af = pd.merge(
    combined_raw,
    snp_af_data_combined[['unique_snp_key', 'AF_Gap']], # Use the combined AF_Gap data
    on='unique_snp_key',
    how='left'
)
combined_with_af['AF_Gap'] = combined_with_af['AF_Gap'].fillna(0) # Fill NaN AF_Gap with 0 (should be minimal now)

# --- Step 3: Group by unique SNP + GO-BUSCO term and collect all initial_sources & relevant SNP data ---
snp_go_data = combined_with_af.groupby([
    "unique_snp_key",
    "go_busco"
]).agg(
    combined_source=('initial_source', lambda x: ', '.join(sorted(x.unique()))),
    AF_Gap=('AF_Gap', 'first')
).reset_index()

# --- Step 4: Now group by the GO-BUSCO term and the NEW combined_source, calculating confidence metrics ---
summary_by_combined_source = (
    snp_go_data.groupby(["go_busco", "combined_source"])
    .agg(
        n_snps=('unique_snp_key', 'count'), # Count unique SNPs
        sum_AF_Gap=('AF_Gap', 'sum'),       # Sum of AF_Gaps for all SNPs in this group
        mean_AF_Gap=('AF_Gap', 'mean'),     # Mean AF_Gap for SNPs in this group
        max_AF_Gap=('AF_Gap', 'max')        # Max AF_Gap among SNPs in this group
    )
    .reset_index()
    .sort_values(["sum_AF_Gap", "n_snps"], ascending=[False, False])
)

print("\n--- GO Term / Gene Summary with Confidence Metrics (Top 20) ---")
print(summary_by_combined_source.head(20))

# --- Step 5: Save the new summary table ---
summary_by_combined_source.to_csv(OUTPUT_SUMMARY_FILE, sep="\t", index=False)
print(f"\nComprehensive summary saved to: {OUTPUT_SUMMARY_FILE}")

print("\nAnalysis complete.")
