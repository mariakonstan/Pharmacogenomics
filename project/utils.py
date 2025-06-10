import os
import time
import pysam
import logging
import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm

from rapidfuzz import fuzz
from Bio import Entrez, SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .settings import EMAIL, FASTA_DIR, FASTA_URLS

Entrez.email = EMAIL
logger = logging.getLogger(__name__)


def hamming_identity(s1, s2):
    """
    Returns percent identity between two equal-length strings.
    """
    if len(s1) != len(s2):
        return 0.0
    matches = sum(c1 == c2 for c1, c2 in zip(s1, s2))
    return matches / len(s1)


def fuzz_similarity(s1, s2):
    """
    Returns rapidfuzz ratio between two non equal-length strings.
    """
    return fuzz.ratio(s1, s2) / 100
    

def clean_sequence(seq: str):
    """
    Removes line-breaks from sequence.
    """
    return seq.replace("\n", "")


def get_gene_sequence(accession, start=0, end=0, strand="", db="nucleotide"):
    """
    Returns Gene sequence from NCBI.
    """
    with Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
        ix_start = start - 1 # Convert to 0-based indexing for Python slicing
        ix_end = end  # Python slices are exclusive at end
        seq = record.seq[ix_start:ix_end]
        return str(seq).upper()


def load_gene_sequence(file_path):
    """
    Returns Gene sequence from local data file.
    """
    with open(file_path, 'r') as f:
        return clean_sequence(f.read()).upper()

    
def download_fasta():
    """
    Downloads necessary input files.
    """
    for url in FASTA_URLS:
        filename = url.split("/")[-1]
        destination = FASTA_DIR / filename
        if not destination.exists():
            logger.info(f"Downloading {url} to {destination}")
            subprocess.run(["wget", "-O", str(destination), url], check=True)
        else:
            logger.info(f"File already exists: {destination}")    
   

def convert_ftp_to_https(url: str) -> str:
    if url.startswith("ftp://"):
        return url.replace("ftp://", "https://", 1)
    return url
   
   
def safe_open_alignment(url, ref_fasta, mode="rc", retries=3, delay=2):
    try:
        return pysam.AlignmentFile(url, mode=mode, reference_filename=ref_fasta)
    except Exception as e:
        if retries > 0:
            logger.warning(f"{url} > connection busy: Retrying... (remaining tries: {retries})")
            time.sleep(delay)
            return safe_open_alignment(url, ref_fasta, mode, retries-1, delay)
        

def plot_total_summary(df, output_dir, folder=None):
    
    cols_skipped = ["Sample", "Population", "Response", "Metabolizer", "Side-Effects"]
    dir = os.path.join(output_dir, folder) if folder else output_dir
    os.makedirs(dir, exist_ok=True)
    
    allele_cols = [col for col in df.columns if "Allele" in col]
    gene_cols = [col for col in df.columns if not "Allele" in col and col not in cols_skipped]
    
    df_melted = df.melt(id_vars=["Sample", "Population"], value_vars=gene_cols, var_name="Gene", value_name="Match Score")
    df_allele = df.set_index("Sample")[allele_cols]
    
    """
    For each gene, visualize the distribution of the total gene match scores.
    Use: Spot genes with generally high or low confidence.
    """
    plt.figure(figsize=(12, 6))
    sns.stripplot(data=df_melted, x="Gene", y="Match Score", jitter=True, size=2)
    plt.xticks(rotation=45)
    plt.title("Gene Match Score Distribution Across Samples")
    plt.tight_layout()
    plt.savefig(f"{dir}/gene_match_distr.png", dpi=300, bbox_inches="tight")  
    plt.close()  
    """
    A heatmap of all allele scores for each sample.
    Use: Spot outliers or patterns in allele-level match strengths.    
    """
    plt.figure(figsize=(18, 10))
    sns.heatmap(df_allele, cmap="viridis", cbar_kws={'label': 'Match Score'})
    plt.title("Allele Match Scores per Sample")
    plt.xlabel("Alleles")
    plt.ylabel("Samples")
    plt.tight_layout()
    plt.savefig(f"{dir}/allele_heatmap.png", dpi=300, bbox_inches="tight")  
    plt.close()  
    """
    Compare how  gene match scores differ across populations.
    Use: Identify population-specific gene score trends.
    """
    grouped = df.groupby("Population")[gene_cols].mean().T
    grouped.plot(kind="bar", figsize=(18, 8), width=0.8)
    plt.title("Mean Gene Match Scores by Population")
    plt.ylabel("Average Match Score")
    plt.xlabel("Genes")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{dir}/gene_match_by_pop.png", dpi=300, bbox_inches="tight")  
    plt.close()  
    """
    Compare how allele match scores differ across populations.
    Use: Identify population-specific allele score trends.
    """
    grouped = df.groupby("Population")[allele_cols].mean().T
    grouped.plot(kind="bar", figsize=(18, 8), width=0.8)
    plt.title("Mean Allele Match Scores by Population")
    plt.ylabel("Average Match Score")
    plt.xlabel("Alleles")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{dir}/allele_match_by_pop.png", dpi=300, bbox_inches="tight")  
    plt.close()  
    """
    Reduce dimensionality of all match scores to see sample clustering.
    Use: Detect whether samples group by population.
    """
    features = df.drop(columns=cols_skipped)
    scaled = StandardScaler().fit_transform(features)
    
    pca = PCA(n_components=2)
    components = pca.fit_transform(scaled)
    
    pca_df = pd.DataFrame(components, columns=["PC1", "PC2"])
    pca_df["Population"] = df["Population"].values
    
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Population")
    plt.title("PCA of Sample Genetic Profiles")
    plt.tight_layout()
    plt.savefig(f"{dir}/genetic_profile_PCA.png", dpi=300, bbox_inches="tight")  
    plt.close()  
    """
    Correlation between allele scores across all samples.
    Use: See which alleles tend to co-occur or behave similarly.
    """    
    allele_corr = df[allele_cols].corr()
    plt.figure(figsize=(14, 12))
    sns.heatmap(allele_corr, cmap="coolwarm", center=0, annot=False)
    plt.title("Correlation Matrix of Allele Match Scores")
    plt.tight_layout()
    plt.savefig(f"{dir}/corr_matrix.png", dpi=300, bbox_inches="tight")  
    plt.close()  
    """
    Total match score per sample across all genes.
    Use: Flag samples with generally low match quality (outlier detection).
    """
    df["Total_Match"] = df[gene_cols].sum(axis=1)
    plt.figure(figsize=(10, 5))
    sns.histplot(df["Total_Match"], bins=30)
    plt.title("Distribution of Total Match Scores per Sample")
    plt.xlabel("Total Match Score")
    plt.ylabel("Sample Count")
    plt.tight_layout()
    plt.savefig(f"{dir}/score_sums_distr.png", dpi=300, bbox_inches="tight")  
    plt.close()  

    logger.info(f"Plots saved to: {dir}/")


def plot_per_gene_summary(df, output_dir, tag):
    """
    Specific Allele set plots
    """
    dir = os.path.join(output_dir, tag)
    os.makedirs(dir, exist_ok=True)

    # Detect allele columns
    allele_cols = [col for col in df.columns if "Allele" in col]
    
    # 1. Allele Frequency Distribution
    plt.figure(figsize=(10, 5))
    df[allele_cols].mean().plot(kind='bar', color='steelblue')
    plt.title(f"{tag} - Allele Frequency Distribution")
    plt.ylabel("Mean Match Score")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"{dir}/allele_freq.png")
    plt.close()

    # 2. Genotype Distribution
    plt.figure(figsize=(8, 5))
    df['Genotype'].value_counts().plot(kind='bar', color='orchid')
    plt.title(f"{tag} - Genotype Distribution")
    plt.ylabel("Sample Count")
    plt.xlabel("Genotype")
    plt.tight_layout()
    plt.savefig(f"{dir}/genotype_distribution.png")
    plt.close()

    # 3. Combined Phenotype Bar Plot
    phenotypes = {
        "Metabolizer": df["Metabolizer"].value_counts(),
        "Response": df["Response"].value_counts(),
        "Side-Effects": df["Side-Effects"].map({1: "Yes"}).dropna().value_counts()
    }
    phenotype_df = pd.DataFrame(phenotypes).fillna(0).astype(int)

    phenotype_df.plot(kind='bar', figsize=(10, 6), color=["#4daf4a", "#377eb8", "#e41a1c"])
    plt.title(f"{tag} - Phenotype Summary")
    plt.ylabel("Sample Count")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f"{dir}/phenotype_summary.png")
    plt.close()

    # 4. Allele Score Heatmap
    plt.figure(figsize=(12, min(0.3 * len(df), 10)))
    sns.heatmap(df[allele_cols], cmap="viridis", yticklabels=False)
    plt.title(f"{tag} - Allele Match Score Heatmap")
    plt.tight_layout()
    plt.savefig(f"{dir}/heatmap.png")
    plt.close()

    logger.info(f"Plots saved to: {dir}/")

        