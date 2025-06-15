import os
import logging
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import cycle
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pycirclize import Circos

from project.models import GeneType

logger = logging.getLogger(__name__)


class Plotter:
    def __init__(self, input_csv, output_dir, pops_json=None):
        self.input_csv = input_csv
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.epsilon = 1e-3
        self.pop_groups = None
        if pops_json:
            import json
            with open(pops_json, 'r', encoding='utf-8') as f:
                self.pop_groups = json.load(f)
            
        self.df = pd.read_csv(os.path.join(self.input_csv))


    def save_plot(self, filename):
        filepath = f"{self.output_dir}/{filename}.png"
        plt.savefig(filepath, dpi=300, bbox_inches="tight")  
        logger.info(f"Plot saved: {filepath}")


class SumPlotter(Plotter):
    def __init__(self, input_csv, output_dir, pops_json=None, non_score_cols=[]):
        super().__init__(input_csv, output_dir, pops_json)
        if non_score_cols:
            self.non_score_cols = non_score_cols
        else:
            self.non_score_cols = ["Sample", "Population", "Response", "Metabolizer", "Side-Effects"]
        self.allele_cols = [col for col in self.df.columns if "Allele" in col]
        self.gene_cols = [col for col in self.df.columns if not "Allele" in col and col not in self.non_score_cols]


    def plot_allele_heatmap_by_sample(self, filename=None):
        """
        A heatmap of all allele scores for each sample.
        Use: Spot sample outliers or patterns in allele-level match strengths.    
        """
        df_copy = self.df.copy()
        df_allele = self.df[self.allele_cols].copy()
        df_allele = df_copy.set_index("Sample")[self.allele_cols]
        plt.figure(figsize=(18, 10))
        sns.heatmap(df_allele, cmap="viridis", cbar_kws={'label': 'Match Score'})
        plt.title("Allele Match Scores per Sample")
        plt.xlabel("Alleles")
        plt.ylabel("Samples")
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
    
    
    def plot_allele_heatmap_by_population(self, filename=None):
        """
        A heatmap of all allele score means by population group.
        Use: Spot population outliers or patterns in allele-level match strengths.    
        """
        # Compute mean allele scores by population
        df_allele = self.df[self.allele_cols].copy()
        df_allele["Population"] = self.df["Population"]
        pop_allele_means = df_allele.groupby("Population").mean()
        # Plot heatmap by population
        plt.figure(figsize=(18, 8))
        sns.heatmap(pop_allele_means, cmap="viridis", cbar_kws={'label': 'Mean Match Score'})
        plt.title("Mean Allele Match Scores by Population")
        plt.xlabel("Alleles")
        plt.ylabel("Population")
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
       
    
    def plot_allele_matches_by_population(self, filename=None):
        """
        Compare how allele match scores differ across populations.
        Use: Identify population-specific allele score trends.
        """
        grouped = self.df.groupby("Population")[self.allele_cols].mean().T
        grouped.plot(kind="bar", figsize=(18, 8), width=0.8)
        plt.title("Mean Allele Match Scores by Population")
        plt.ylabel("Average Match Score")
        plt.xlabel("Alleles")
        plt.xticks(rotation=90)
        plt.legend(title="Population", bbox_to_anchor=(0.5, -0.4), loc="upper center", ncol=13)
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
       
    
    def plot_pca_by_population(self, filename=None):
        """
        Reduce dimensionality of all match scores to see sample clustering.
        Use: Detect whether samples group by population.
        """
        features = self.df.drop(columns=self.non_score_cols)
        scaled = StandardScaler().fit_transform(features)
        pca = PCA(n_components=2)
        components = pca.fit_transform(scaled)
        
        pca_df = pd.DataFrame(components, columns=["PC1", "PC2"])
        pca_df["Population"] = self.df["Population"].values
        
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Population", s=50)
        plt.title("PCA of Sample Genetic Profiles")
        plt.legend(title="Population", bbox_to_anchor=(0.5, -0.22), loc="upper center", ncol=10)
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
        
    
    def plot_pca_by_continent(self, filename=None):
        """
        Reduce dimensionality of all match scores to see sample clustering.
        Use: Detect whether samples group by a wider population group (continental).
        """
        assert(self.pop_groups)
        df_copy = self.df.copy()
        df_copy["Continent"] = df_copy["Population"].map(self.pop_groups)
        df_filtered = df_copy.dropna(subset=["Continent"]).copy()
        
        features = df_filtered.drop(columns=self.non_score_cols + ["Continent"])
        features = features.select_dtypes(include=["number"])
    
        scaled = StandardScaler().fit_transform(features)
        pca = PCA(n_components=2)
        components = pca.fit_transform(scaled)
        
        pca_df = pd.DataFrame(components, columns=["PC1", "PC2"])
        pca_df["Continent"] = df_filtered["Continent"].values
        
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Continent", s=50)
        plt.title("PCA of Sample Genetic Profiles (by Continent)")
        plt.legend(title="Population Group", bbox_to_anchor=(0.5, -0.22), loc="upper center", ncol=10)
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
    
        
    def plot_correlation_matrix(self, filename=None):
        """
        Correlation between allele scores across all samples.
        Use: See which alleles tend to co-occur or behave similarly.
        """    
        allele_corr = self.df[self.allele_cols].corr()
        plt.figure(figsize=(14, 12))
        sns.heatmap(allele_corr, cmap="coolwarm", center=0, annot=False)
        plt.title("Correlation Matrix of Allele Match Scores")
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
       
    
    def plot_match_sum_across_genes(self, filename=None):
        """
        Total match score per sample across all genes.
        Use: Flag samples with generally low match quality (outlier detection).
        """
        df_copy = self.df.copy()
        df_copy["Total_Match"] = df_copy[self.gene_cols].sum(axis=1)
        df_adjusted = df_copy["Total_Match"] + self.epsilon
        
        plt.figure(figsize=(10, 5))
        ax = sns.histplot(df_adjusted, bins=10)
        plt.yscale("log")
        plt.grid(True, which="both", axis="y", linestyle="--", linewidth=0.7, alpha=0.7)
        plt.title("Distribution of Total Match Scores per Sample")
        plt.xlabel("Total Match Score")
        plt.ylabel("Sample Count (log)")
        # Add counts above each bin
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.text(
                    p.get_x() + p.get_width() / 2,  # x position: center of the bin
                    height * 1.1,                   # y position: slightly above the bar
                    f"{int(height)}",               # text: count as integer
                    ha='center', va='bottom', fontsize=8
                )    
        
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
        
        
    def plot_response_by_population(self, filename=None):
        # Response distribution 
        plt.figure(figsize=(8, 6))
        sns.countplot(data=self.df, x="Population", hue="Response")
        plt.title("Response Types by Population")
        plt.xticks(rotation=45)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
    
    
    def plot_metabolizer_by_population(self, filename=None):
        # Metabolizer distribution
        plt.figure(figsize=(8, 6))
        sns.countplot(data=self.df, x="Population", hue="Metabolizer")
        plt.title("Metabolizer Types by Population")
        plt.xticks(rotation=45)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
    
    
    def plot_response_ratio_by_population(self, filename=None):
        data = self.df.groupby(['Population', "Response"]).size().unstack(fill_value=0)
        data_percent = data.div(data.sum(axis=1), axis=0)
        data_percent.plot(kind='bar', stacked=True, colormap='Set2', figsize=(10,6))
        plt.ylabel('Proportion')
        plt.title('Response by Population')
        plt.legend(title="Response")
        plt.xticks(rotation=45)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
    

class GeneGroupPlotter(Plotter):
    def __init__(self, group, input_csv, output_dir, pops_json=None):
        super().__init__(input_csv, output_dir, pops_json)
        self.group = group
        self.allele_cols = [col for col in self.df.columns if "Allele" in col]


    def plot_heatmap(self, filename=None):
        # Allele Score Heatmap
        plt.figure(figsize=(12, min(0.3 * len(self.df), 10)))
        sns.heatmap(self.df[self.allele_cols], cmap="viridis", yticklabels=False)
        plt.title(f"{self.group} - Allele Match Score Heatmap")
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
   
    
    def plot_frequency_distribution(self, filename=None):
        # Allele Frequency Distribution
        plt.figure(figsize=(10, 5))
        self.df[self.allele_cols].mean().plot(kind='bar', color='steelblue')
        plt.title(f"{self.group} - Allele Frequency Distribution")
        plt.ylabel("Mean Match Score")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()

    
    def plot_genotype_distribution(self, filename=None):
        # Genotype Distribution
        plt.figure(figsize=(8, 5))
        self.df['Genotype'].value_counts().plot(kind='bar', color='orchid')
        plt.title(f"{self.group} - Genotype Distribution")
        plt.ylabel("Sample Count")
        plt.xlabel("Genotype")
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()
       
        
    def plot_phenotype_distribution(self, filename=None):
        # Combined Phenotype Bar Plot
        phenotypes = {
            "Metabolizer": self.df["Metabolizer"].value_counts(),
            "Response": self.df["Response"].value_counts(),
            "Side-Effects": self.df["Side-Effects"].map({1: "Yes"}).dropna().value_counts()
        }
        phenotype_df = pd.DataFrame(phenotypes).fillna(0).astype(int)
        phenotype_df.plot(kind='bar', figsize=(10, 6), color=["#4daf4a", "#377eb8", "#e41a1c"])
        plt.title(f"{self.group} - Phenotype Summary")
        plt.ylabel("Sample Count")
        plt.xticks(rotation=0)
        plt.tight_layout()
        if filename:
            self.save_plot(filename)
            plt.close()  
        else:
            plt.show()


def plot_circos_genome(gene_data, matches, output_dir=None, folder=None):
    dir = os.path.join(output_dir, folder) if folder else output_dir
    os.makedirs(dir, exist_ok=True)
    
    # Step 1: Build chromosome size map
    chrom_sizes = defaultdict(int)
    for g in gene_data:
        if not g.chrom:
            continue
        chrom = g.chrom.replace("chr", "")
        end = g.end if g.type.name == "GENE" else g.pos
        chrom_sizes[chrom] = max(chrom_sizes[chrom], end + 1_000_000)

    # Step 2: Initialize Circos plot
    circos = Circos(dict(chrom_sizes), space=2)

    # Track handle for each chromosome
    tracks = {}
    for chrom in chrom_sizes:
        sector = next(s for s in circos.sectors if s.name == chrom)
        tracks[chrom] = sector.add_track((60, 80))  # Outer ring

    # Default tag colors if none provided
    tags = list(set(g.tag for g in gene_data if g.tag))
    palette = cycle(sns.color_palette("tab10", len(tags)).as_hex())
    tag_colors = {tag: next(palette) for tag in tags}

    # Step 3: Draw genes/alleles in track
    for g in gene_data:
        chrom = g.chrom.replace("chr", "")
        tr = tracks[chrom]

        start = g.start if g.type == GeneType.GENE else g.pos
        end = g.end if g.type == GeneType.GENE else g.pos + 1
        color = tag_colors.get(g.tag, "lightgrey")

        tr.axis(fc=color, ec="black", alpha=0.6)
        tr.rect(start, end, fc=color, ec="black", lw=0.3)
        tr.text(g.name, (start + end) / 2, r=75, size=4, ha="center", va="center")

    # Step 4: Add links from matches (within same chrom)
    for m in matches:
        g = next((g for g in gene_data if g.name == m.gene_name), None)
        if not g or not g.chrom:
            continue

        chrom = g.chrom.replace("chr", "")
        pos1 = g.start if g.type.name == "GENE" else g.pos
        pos2 = g.end if g.type.name == "GENE" else g.pos + 1

        score = m.score if m.score is not None else 0.5
        color = tag_colors.get(g.tag, "grey")
        alpha = 0.3 + 0.7 * score
        linewidth = 0.5 + score * 2

        circos.link((chrom, pos1, pos2), (chrom, pos1, pos2), color=color, lw=linewidth, alpha=alpha)

    # Step 5: Render plot
    fig = circos.plotfig()
    fig.savefig(f"{dir}/circos_genome.png", dpi=300, bbox_inches="tight")

