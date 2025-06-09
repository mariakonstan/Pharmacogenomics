import os
import copy
import json
import logging
import pandas as pd
import concurrent.futures as ftr

from enum import Enum
from collections import Counter

import project.utils as utils
from project.settings import FUZZY_IDENTITY_THRESHOLD

logger = logging.getLogger(__name__)


class SampleData:
    def __init__(self, name, url, pop=""):
        self.name = name
        self.url = url
        self.pop = pop


class Match:
    def __init__(self, sample_id, sample_pop, gene_name, match_type, score=None):
        """
        Match holds sample matching type and score.
        Args:
            sample_id (str): sample name or identifier
            sample_population (str): sample population 
            gene_name (str): gene/allele name or identifier
            match_type (str): characterization of matching (e.g. gene, gene_fuzzy, allele)
            score (float): 1.0 for complete match (< 1.0 fuzzy match)
        """
        self.sample_id = sample_id
        self.sample_population = sample_pop
        self.gene_name = gene_name 
        self.match_type = match_type  
        self.score = score  # could be frequency %, match % (max: 1.0)

    def to_dict(self):
        return {
            "sample_id": self.sample_id,
            "sample_population": self.sample_population,
            "gene_name": self.gene_name,
            "match_type": self.match_type,
            "score": self.score
        }


class GeneType(Enum):
    NONE = ""
    GENE = "GENE"
    ALLELE = "ALLELE"


class GeneData:
    def __init__(self, name, chrom="", start=0, end=0):
        """
        GeneData holds genomic region info.
        Args:
            name (str): gene name or identifier
            chrom (str): chromosome name (e.g. 'chr17')
            start (int): 1-based start position
            end (int): end position
        """
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.type = GeneType.NONE
        self.pos = None # if defined it's a gene allele
        self.seq = "" # nucleotide sequence (len should be end-start)
        self.tag = "" # some label to group gene-data
        self.attr = "" # characterization depending on drug response
        self.comments = "" 
        # Per sample evaluation 
        self.match = None # sample match

    @staticmethod
    def create(g:dict, data_dir=""):
        gene = GeneData(name=g["name"], chrom=g["chrom"])
        gene.pos = g.get("pos", "")
        gene.tag = g.get("tag", "")
        gene.attr = g.get("attr", "")
        gene.comments = g.get("comments", "")
        seq_accession = g.get("accession", "")
        seq_file = g.get("file", "")
        
        if gene.pos:
            gene.type = GeneType.ALLELE 
            gene.start=gene.pos + g.get("start", 0) # start offset from pos
            gene.end=gene.pos + g.get("end", 0) # end offset from pos
            gene.seq = g.get("base", "")
        else:
            gene.type = GeneType.GENE
            gene.start = g.get("start", 0) # defined start pos
            gene.end = g.get("end", 0) # defined end pos
            if seq_file: # File paths already saved from Entrez 
                file_path = os.path.join(data_dir, seq_file)
                if os.path.exists(file_path):
                    gene.seq = utils.load_gene_sequence(file_path)
                else:
                    if seq_accession:
                        gene.seq = utils.get_gene_sequence(seq_accession, gene.start, gene.end)
                        with open(file_path, "w") as file:
                            logger.info(f"Saving sequence locally for: {gene.name}")
                            file.write(gene.seq)
                    else:
                        logger.warning(f"No sequence accession number given for gene: {gene.name}")
            elif seq_accession:
                gene.seq = utils.get_gene_sequence(seq_accession, gene.start, gene.end)
            
        return gene

    def query_sample(self, sample: SampleData, ref_fasta):
        """
        Query a remote CRAM file for reads overlapping this gene region.
        Args:
            sample (SampleData): sample related info (Name, CRAM file URL, etc.)
            ref_fasta (str): local reference fasta path
        """
        self.match = None # Reset match
        try:
            bamfile = utils.safe_open_alignment(sample.url, ref_fasta, mode='rc', retries=2, delay=2)
            logger.info(f"{sample.name}: matching '{self.name}' sequence (size:{len(self.seq)})")
            # Prepare allele check (if self.pos is defined)
            pileup_checked = False
            if self.pos:
                for pileupcolumn in bamfile.pileup(self.chrom, self.pos - 1, self.pos, truncate=True):
                    if pileupcolumn.pos == self.pos - 1:
                        bases = [
                            pileupread.alignment.query_sequence[pileupread.query_position].upper()
                            for pileupread in pileupcolumn.pileups
                            if pileupread.query_position is not None
                        ]
                        pileup_checked = True
                        base_counts = Counter(bases)
                        total = sum(base_counts.values())
            
                        expected_base = self.seq.upper()
                        allele_count = base_counts.get(expected_base, 0)
                        allele_freq = allele_count / total if total else 0.0
            
                        if allele_count > 0:
                            self.match = Match(
                                sample_id=sample.name, 
                                sample_pop=sample.pop,
                                gene_name=self.name,
                                match_type="allele", 
                                score=round(allele_freq, 4)
                            )
                        break  # Done checking allele at this position                
                if not pileup_checked:
                    logger.warning(f"No pileup column found at {self.chrom}:{self.pos} for {self.name}")
            else:
                # Iterate once through reads for sequence matching
                ix_start = self.start-1
                ix_end = self.end
                consensus = [] # Reconstruct consensus sequence from pileup
                for column in bamfile.pileup(self.chrom, ix_start, ix_end, truncate=True):
                    if ix_start <= column.reference_pos < ix_end:
                        bases = [
                            pileupread.alignment.query_sequence[pileupread.query_position]
                            for pileupread in column.pileups
                            if pileupread.query_position is not None and not pileupread.is_del
                        ]
                        if bases:
                            most_common_base, _ = Counter(bases).most_common(1)[0]
                            consensus.append(most_common_base)
                        else:
                            consensus.append('-')  # No coverage
                
                reconstructed_seq = ''.join(consensus).upper()
                
                if len(reconstructed_seq) == len(self.seq):
                    identity = utils.hamming_identity(reconstructed_seq, self.seq)
                    if identity == 1.0:
                        self.match = Match(
                            sample_id=sample.name, 
                            sample_pop=sample.pop,
                            gene_name=self.name,
                            match_type="gene", 
                            score=1.0
                        )
                    elif identity >= FUZZY_IDENTITY_THRESHOLD:
                        self.match = Match(
                            sample_id=sample.name, 
                            sample_pop=sample.pop,
                            gene_name=self.name,
                            match_type="gene_fuzzy", 
                            score=round(identity, 4)
                        )
                else:
                    logger.warning(f"{self.name} at sample: {sample.name}: sequence lengths don't match (query: {len(self.seq)}, read:{len(reconstructed_seq)}). Using rapidfuzz ratio.")
                    similarity = utils.fuzz_similarity(reconstructed_seq, self.seq)
                    if similarity >= FUZZY_IDENTITY_THRESHOLD:
                        self.match = Match(
                            sample_id=sample.name, 
                            sample_pop=sample.pop,
                            gene_name=self.name,
                            match_type="gene_fuzzy_indel", 
                            score=round(similarity, 4)
                        )
            bamfile.close()                
    
        except Exception as e:
            logger.error(f"Error querying {sample.url} for {self.name} in sample {sample.name}: {e}")

    def summarize(self):
        """
        Summarize query results 
        """
        return {
            "gene": self.name,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "tag": self.tag,
            "attr": self.attr,
            "comments": self.comments,
            "match": self.match.to_dict() if self.match else {}
        }
        

class GeneDataGroup:
    def __init__(self, tag, gene:GeneData, alleles=[]):
        self.tag = tag
        self.gene = gene
        self.alleles = alleles


class SampleScanner:
    def __init__(self):
        self.index_file = None
        self.ref_fasta_file = None
        self.output_dir = None
        self.output_dir_samples = None
        self.output_dir_figures = None
        self.genes_json = None
        self.gene_data = []
        self.gene_groups = {}
        self.limit = 0 # sample scan limit (for debug)
        self.start = 0 # index to continue from (on scan broken)

    @staticmethod    
    def read_sample_matches(file_path):
        sample_matches = {}
        with open(file_path, 'r', encoding='utf-8') as f:
            sample_data = json.load(f)
            for entry in sample_data:
                match = None
                key = entry.get("gene", "")
                entry_match = entry.get("match", {})
                if len(entry_match):
                    match = Match(
                        sample_id=entry_match.get("sample_id", ""),
                        sample_pop=entry_match.get("sample_population", ""),
                        gene_name=entry_match.get("gene_name", ""),
                        match_type=entry_match.get("match_type", ""),
                        score=entry_match.get("score", 0)
                    )
                sample_matches[key] = match
        return sample_matches
    
    @staticmethod    
    def read_sample_matches_lite(file_path):
        sample_matches = {}
        with open(file_path, 'r', encoding='utf-8') as f:
            match_data = json.load(f)
            for entry in match_data:
                match = Match(
                    sample_id=entry.get("sample_id", ""),
                    sample_pop=entry.get("sample_population", ""),
                    gene_name=entry.get("gene_name", ""),
                    match_type=entry.get("match_type", ""),
                    score=entry.get("score", 0)
                )
                sample_matches[match.gene_name] = match
        return sample_matches
    
    def add_gene_data(self, gene: GeneData):
        n_sequence = len(gene.seq)
        if n_sequence:
            self.gene_data.append(gene)
            n_region = gene.end - gene.start + 1
            if n_region and n_sequence > n_region:
                logger.warning(f"Gene {gene.name} seq length ({n_sequence}) is larger than region length ({n_region})")
            logger.info(f"Gene data loaded: '{gene.name}'")
        else:
            logger.warning(f"Sequence data missing for '{gene.name}'")
        
    def load_gene_data(self):
        data_dir = os.path.dirname(os.path.abspath(self.genes_json))
        with open(self.genes_json) as f:
            gene_list = json.load(f)
        self.gene_data = []
        self.gene_groups = {}
        for g in gene_list:
            gene = GeneData.create(g, data_dir)
            self.add_gene_data(gene)
            alleles = g.get("alleles", [])
            tags = {}
            for a in alleles:
                allele = GeneData.create(a, data_dir)
                self.add_gene_data(allele)
                if allele.tag not in tags.keys():
                    tags[allele.tag] = []
                tags[allele.tag].append(allele)
            # Group alleles (only) by tag
            for tag, alleles in tags.items():
                self.gene_groups[tag] = GeneDataGroup(tag, gene, alleles)
            
    def process_sample(self, sample_row):
        """
        Process one sample (row) and return a dict of sample results (score per gene/allele)
        """
        sample_name = sample_row["SAMPLE_NAME"]
        sample_pop = sample_row["POPULATION"]
        sample_url = utils.convert_ftp_to_https(sample_row["ENA_FILE_PATH"])
        sample_result = {"sample": sample_name}
        sample = SampleData(name=sample_name, url=sample_url, pop=sample_pop)
        sample_summaries = []

        # Make deep copies of gene_data to avoid shared state during parallel runs
        sample_genes = [copy.deepcopy(gene) for gene in self.gene_data]
    
        try:
            for gene in sample_genes:
                gene.query_sample(sample, self.ref_fasta_file)
                summary = gene.summarize()
                sample_summaries.append(summary)
                """ TODO: export gene-match dict only instead of gene-summary
                if gene.match:
                    sample_summaries.append(gene.match.to_dict())
                """
                sample_result[gene.name] = gene.match.score if gene.match else 0
        except Exception as e:
            logger.error(f"Error processing {sample_name}: {e}")
                
        # Delete temp CRAI files
        try:
            crai_path = sample.url.split("/")[-1] + ".crai"  
            if os.path.exists(crai_path):
                os.remove(crai_path)
        except Exception as e:
            logger.warning(f"Could not delete CRAI index for {sample.name}: {e}")                
        
        self.export_sample_data(sample_summaries, sample_name)
        return sample_result

    def export_csv(self, csv_data, output_csv=None):
        if output_csv is None:
            logger.warning(f"No file name given to export CSV data!")
            return
        path = os.path.join(self.output_dir, f"{output_csv}.csv")
        flat_df = pd.DataFrame(csv_data)
        flat_df.to_csv(path, index=False)
        logger.info(f"CSV summary saved to {path}")
    
    def export_sample_data(self, json_data, sample_name):
        if sample_name is None:
            return
        path = os.path.join(self.output_dir_samples, f"{sample_name}.json")
        with open(path, "w") as f:
            json.dump(json_data, f, indent=2)
        logger.info(f"Sample summary saved to: {path}")

    def _setup(self):
        assert(self.index_file)
        assert(self.output_dir)
        assert(self.ref_fasta_file)
        assert(len(self.gene_data))

    def run(self):
        self._setup()
        df = pd.read_csv(self.index_file, sep="\t")
        n = len(df)

        if self.start:
            df = df.iloc[self.start-1:]
        if self.limit:
            df = df.iloc[:self.limit]

        for idx, row in df.iterrows():
            logger.info(f"Processing sample {row['SAMPLE_NAME']} ({idx + 1}/{n})")
            self.process_sample(row)

    def run_parallel(self, max_workers=4):
        self._setup()
        df = pd.read_csv(self.index_file, sep="\t")
        n = len(df)
        
        if self.start:
            df = df.iloc[self.start-1:]
        if self.limit:
            df = df.iloc[:self.limit]

        logger.info(f"Running in parallel with {max_workers} workers")

        with ftr.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Map each sample row to process_sample, returns list of futures
            futures = [executor.submit(self.process_sample, row) for _, row in df.iterrows()]

            for i, future in enumerate(ftr.as_completed(futures)):
                try:
                    sample_result = future.result()
                    logger.info(f"Completed sample {sample_result['sample']} ({i+1}/{n})")
                except Exception as e:
                    logger.error(f"Error processing sample in parallel: {e}")
    
    def export_summaries(self):
        total_summary = []
        group_summaries = {}
        for tag, group in self.gene_groups.items():
            group_summaries[tag] = []
            
        df = pd.read_csv(self.index_file, sep="\t")
        
        folder_path = self.output_dir_samples 
        for filename in os.listdir(folder_path):
            if filename.endswith('.json'):
                sample_name = os.path.splitext(filename)[0]
                pop = df.loc[df["SAMPLE_NAME"] == sample_name, "POPULATION"].values
                population = pop[0] if pop.size > 0 else ""
                
                file_path = os.path.join(folder_path, filename)
                sample_matches = SampleScanner.read_sample_matches(file_path)
                
                # add total summary entry
                sample_entry = {"Sample": sample_name, "Population": population}
                for key, match in sample_matches.items():
                    sample_entry[key] = match.score if match else 0
                total_summary.append(sample_entry)
                
                # update gene_data matches for this sample
                for gene in self.gene_data:
                    gene.match = sample_matches.get(gene.name, None)
                    if gene.type == GeneType.GENE and gene.match is None:
                        logger.warning(f"{sample_name} - Gene: {gene.name}: Gene match not found!")
                        
                # add group-summaries entry
                for tag, group in self.gene_groups.items():
                    group_entry = {"Sample": sample_name, "Population": population}
                    alleles_found = {}
                    genotype = ""
                    response = ""
                    metabolizer = ""
                    side_effects = 0
                
                    # Collect matching allele scores and attributes
                    for allele in group.alleles:
                        score = allele.match.score if allele.match else 0
                        group_entry[allele.name] = score
                
                        if score > 0.3:
                            alleles_found[allele.seq] = allele
                            if allele.attr == "XX":
                                side_effects = 1
                
                    num_alleles = len(alleles_found)
                
                    # Handle cases based on how many alleles passed the score threshold
                    if num_alleles == 0:
                        logger.warning(f"{sample_name} - Gene: {group.gene.name} ({tag}): Allele match not found!")
                    elif num_alleles == 1:
                        seq, allele = next(iter(alleles_found.items()))
                        genotype = seq * 2
                        if allele.attr.endswith("M"):
                            metabolizer = allele.attr
                        elif allele.attr.endswith("R"):
                            response = allele.attr
                    elif num_alleles == 2:
                        # Sort for consistent genotype ordering
                        items = list(alleles_found.items())
                        seq_1, allele_1 = items[0]
                        seq_2, allele_2 = items[-1]
                        wt_allele = None
                        gt = []
                        if "wild type" in allele_1.comments:
                            wt_allele = allele_1
                            gt = [seq_1, seq_2]
                        elif "wild type" in allele_2.comments:
                            wt_allele = allele_2
                            gt = [seq_2, seq_1]
                        else:
                            gt = sorted([seq_1, seq_2])
                        genotype = "".join(gt)
                        
                        # Determine dominant attribute (prefer non-N or wt)
                        if not allele_1.attr.startswith("N"):
                            dom_allele = allele_1
                        elif not allele_2.attr.startswith("N"):
                            dom_allele = allele_2
                        else:
                            dom_allele = wt_allele if wt_allele else allele_1
                        
                        if dom_allele.attr.endswith("M"):
                            metabolizer = dom_allele.attr
                        elif dom_allele.attr.endswith("R"):
                            response = dom_allele.attr
                    else:
                        logger.warning(f"{sample_name} - Gene: {group.gene.name} ({tag}): More than 2 alleles found!")

                    group_entry["Genotype"] = genotype 
                    group_entry["Response"] = response 
                    group_entry["Metabolizer"] = metabolizer 
                    group_entry["Side-Effects"] = side_effects 
                    group_summaries[tag].append(group_entry)
        
        for tag, group_summary in group_summaries.items():
            self.export_csv(group_summary, f"summary_{tag}")
        self.export_csv(total_summary, "summary_total")
    
    def plot_figures(self):
        fig_dir = self.output_dir_figures
        df = pd.read_csv(os.path.join(self.output_dir, "summary_total.csv"))
        utils.plot_gene_match_distribution(df, os.path.join(fig_dir, "fig_gene_distr.png"))
        utils.plot_gene_match_by_population(df, os.path.join(fig_dir, "fig_gene_pop.png"))
        utils.plot_allele_match_by_population(df, os.path.join(fig_dir, "fig_allele_pop.png"))
        utils.plot_allele_match_corr_matrix(df, os.path.join(fig_dir, "fig_allele_corr.png"))
        utils.plot_allele_heatmap(df, os.path.join(fig_dir, "fig_allele_heat.png"))
        utils.plot_genetic_profile_PCA(df, os.path.join(fig_dir, "fig_genetic_profile.png"))
        utils.plot_total_match_distribution(df, os.path.join(fig_dir, "fig_total_distr.png"))