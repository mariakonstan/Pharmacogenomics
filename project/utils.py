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
        

def plot_gene_match_distribution(df, save_path=None):
    """
    For each gene, visualize the distribution of the total gene match scores.
    Use: Spot genes with generally high or low confidence.
    """
    genes = [col for col in df.columns if not "Allele" in col and col not in ["Sample", "Population"]]
    df_melted = df.melt(id_vars=["Sample", "Population"], value_vars=genes, var_name="Gene", value_name="Match Score")
    plt.figure(figsize=(12, 6))
    sns.stripplot(data=df_melted, x="Gene", y="Match Score", jitter=True, size=2)
    plt.xticks(rotation=45)
    plt.title("Gene Match Score Distribution Across Samples")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")  
        plt.close()  
    else:
        plt.show()


def plot_allele_heatmap(df, save_path=None):
    """
    A heatmap of all allele scores for each sample.
    Use: Spot outliers or patterns in allele-level match strengths.    
    """
    allele_cols = [col for col in df.columns if "Allele" in col]
    allele_df = df.set_index("Sample")[allele_cols]
    plt.figure(figsize=(18, 10))
    sns.heatmap(allele_df, cmap="viridis", cbar_kws={'label': 'Match Score'})
    plt.title("Allele Match Scores per Sample")
    plt.xlabel("Alleles")
    plt.ylabel("Samples")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")  
        plt.close()  
    else:
        plt.show()


def plot_gene_match_by_population(df, save_path=None):
    """
    Mean gene scores across populations.
    Use: Identify if certain genes are less well-matched in certain populations.
    """
    genes = [col for col in df.columns if not "Allele" in col and col not in ["Sample", "Population"]]
    df_melted = df.melt(id_vars=["Sample", "Population"], value_vars=genes, var_name="Gene", value_name="Match Score")
    plt.figure(figsize=(12, 6))
    sns.swarmplot(data=df_melted, x='Gene', y='Match Score', hue="Population", size=2)
    plt.title("Gene Match Scores by Population")
    plt.xticks(rotation=45)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")  
        plt.close()  
    else:
        plt.show()


def plot_allele_match_by_population(df, save_path=None):
    """
    Compare how allele match scores differ across populations.
    Use: Identify population-specific allele score trends.
    """
    allele_cols = [col for col in df.columns if "Allele" in col]
    grouped = df.groupby("Population")[allele_cols].mean().T
    grouped.plot(kind="bar", figsize=(18, 8), width=0.8)
    plt.title("Mean Allele Match Scores by Population")
    plt.ylabel("Average Match Score")
    plt.xlabel("Alleles")
    plt.xticks(rotation=90)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")  
        plt.close()  
    else:
        plt.show()


def plot_genetic_profile_PCA(df, save_path=None):
    """
    Reduce dimensionality of all match scores to see sample clustering.
    Use: Detect whether samples group by population.
    """
    features = df.drop(columns=["Sample", "Population"])
    scaled = StandardScaler().fit_transform(features)
    
    pca = PCA(n_components=2)
    components = pca.fit_transform(scaled)
    
    pca_df = pd.DataFrame(components, columns=["PC1", "PC2"])
    pca_df["Population"] = df["Population"].values
    
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Population")
    plt.title("PCA of Sample Genetic Profiles")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")  
        plt.close()  
    else:
        plt.show()


def plot_allele_match_corr_matrix(df, save_path=None):
    """
    Correlation between allele scores across all samples.
    Use: See which alleles tend to co-occur or behave similarly.
    """    
    allele_cols = [col for col in df.columns if "Allele" in col]
    allele_corr = df[allele_cols].corr()
    plt.figure(figsize=(14, 12))
    sns.heatmap(allele_corr, cmap="coolwarm", center=0, annot=False)
    plt.title("Correlation Matrix of Allele Match Scores")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")  
        plt.close()  
    else:
        plt.show()


def plot_total_match_distribution(df, save_path=None):
    """
    Total match score per sample across all genes.
    Use: Flag samples with generally low match quality (outlier detection).
    """
    genes = [col for col in df.columns if not "Allele" in col and col not in ["Sample", "Population"]]
    df["Total_Match"] = df[genes].sum(axis=1)
    plt.figure(figsize=(10, 5))
    sns.histplot(df["Total_Match"], bins=30)
    plt.title("Distribution of Total Match Scores per Sample")
    plt.xlabel("Total Match Score")
    plt.ylabel("Sample Count")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")  
        plt.close()  
    else:
        plt.show()
        