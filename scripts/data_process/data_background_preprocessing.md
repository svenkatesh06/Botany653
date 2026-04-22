# DATA BACKGROUND

The primary reference papers for this project are:
1.  *[Genomic basis for the convergent evolution of electric organs, Gallant et al. (2014)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5541775/)*. Gallant et al. studied how a myogenic electric organ, a muscle-derived organ, evolved repeatedly across distantly related fish lineages. The paper assembled the genome of the electric eel (*Electrophorus electricus*) and compared electric organ (EO) and skeletal muscle transcriptomes across five electric-fish species representing multiple independent origins of electrogenesis. A major result of that study was the identification of EO-associated gene sets, including genes in **cluster 9** that were highly upregulated in EO and genes in **cluster 1** that were downregulated in EO relative to muscle.

For this project, I am **not attempting to reproduce the 2014 transcriptome analysis exactly**. Instead, I am using the 2014 paper as the biological starting point, specifically its EO-associated gene sets from cluster 1 and cluster 9, and building an updated comparative sequence analysis around them.

2.  *[Genomic Evidence for Convergent Molecular Adaptation in Electric Fishes, Wang and Yang (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7952227/)*. That study moved from transcriptomic comparisons toward comparative genomics and found evidence for convergent molecular evolution across electric-fish lineages, including accelerated protein evolution, positive selection, and convergent amino acid substitutions. Across 4,657 orthologs, they reported 702 positively selected genes and highlighted functional categories such as cell membrane structure, ion channels, and transmembrane transporters, with **scn4aa** standing out as an especially relevant candidate for electrical and muscular activity.

Together, these 2 papers suggest a useful updated project question:

> **How do modern orthologous coding sequences of EO-associated genes compare across a selected set of electric fish genomes?**

This reframes the original transcriptome-centered question into a genome-based comparative analysis of coding sequences while still staying anchored to EO-associated genes identified by Gallant et al. (2014). In other words:

> **This project is an updated genome-based comparative analysis of EO-associated cluster genes, inspired by the original paper.**

> **I used the gene sets defined in cluster 1 and cluster 9 of the original EO transcriptomic study, but performed an updated comparative sequence analysis using currently available annotated genomes from selected electric fish species.**


# PROJECT DATA AND BIG PICTURE

## Project rationale

The original 2014 species set was ideal for studying EO versus muscle expression at the time. However,some of the original species were analyzed through older transcriptome assemblies, and in several cases gene naming, transcript models, and sequence-to-gene linking are incomplete or inconsistent across species.

Because of that, I am **not using the 2014 species set unchanged**. Instead, I am prioritizing:

1. currently available annotated genomes,
2. species for which CDS/transcript products can be obtained in a relatively uniform format,
3. retention of the broader evolutionary question of electric-organ convergence.

So, rather than reproducing the original RNA-seq dataset directly, I am using the **cluster 1 and cluster 9 gene names** from Gallant et al. (2014) as the biologically meaningful input gene set, and then retrieving modern orthologous coding sequences from currently available genome annotations.

## Big-picture biological question

This project is motivated by two related questions:

1. **Do EO-associated genes identified in the original transcriptomic study show recognizable coding-sequence conservation or divergence across electric-fish genomes?**
2. **Can a modern comparative CDS-based view of these genes reveal patterns consistent with convergent molecular adaptation across electric fish lineages?**

The first question is anchored directly in Gallant et al. (2014), while the second is motivated by Wang and Yang (2021), who showed that electric fishes from multiple lineages exhibit convergent molecular signals at the genome level, including positive selection and convergent substitutions in genes relevant to electrical and muscular function.

### Planned project design

My plan is to use the EO-associated gene names from **cluster 1** and **cluster 9** of Gallant et al. (2014), identify orthologous coding sequences for those genes in a selected set of electric-fish genomes, and then perform per-gene comparative analyses on the CDS sequences. These analyses may include gene-wise alignments and exploratory phylogenetic comparisons across the selected species.

## Species choice

The working logic is:

- include representative electric fishes from more than one lineage,
- include at least one non-electric annotated outgroup,
- avoid species that require heavy manual rescue from older transcriptome-only resources.

With that in mind, the updated species list is:

1. *Electrophorus electricus* - **electric eel**
2. *Brachyhypopomus gauderio* - **brown ghost knifefish**
3. *Paramormyrops kingsleyae* - **weakly electric mormyrid**
4. *Brienomyrus brachyistius* - **baby whale fish**
5. *Ictalurus punctatus* - **channel catfish** *(non-electric outgroup)*

This species list was revised after excluding *Electrophorus voltai* from the project because its current annotation was not convenient for linking sequences to common gene names in a clean and uniform way for downstream analysis. *Brachyhypopomus gauderio* was added instead because it remains biologically relevant as an electric gymnotiform and has a more practical annotation framework for an updated coding-sequence-based comparison.

Here, *Electrophorus electricus* and *Brachyhypopomus gauderio* represent the gymnotiform side of electric fishes, *Paramormyrops kingsleyae* and *Brienomyrus brachyistius* represent the African weakly electric mormyrid side, and *Ictalurus punctatus* provides a non-electric reference point for comparative analysis. This keeps the project focused on electric-fish diversity.

# ENVIRONMENT SETUP

## Install the necessary packages for project in conda environment

I initially tried setting up the conda environment on my Windows system, but some of the required tools are much easier to manage on Linux. I then tried WSL Ubuntu, but it was slow for this workflow, so I finally set up the environment on the WID server.

```bash
conda create -p "<...path_to_env>/e_phylo" \
 python=3.11 \
 sra-tools \
 mafft muscle clustalw \
 r-base r-essentials \
 r-ape r-adegenet r-phangorn \
 bioconductor-ggtree

conda activate envs/e_phylo
```

# DATA : Download and pre-process

## Download annotated CDS FASTA files and annotations

I have obtained the CDS FASTA file and the corresponding annotations for each of 5 species from the NCBI genome. The steps I followed to download the relevant files are menetioned below with reference to 1 species. The same steps were followed for the other 4.

1. Go to https://www.ncbi.nlm.nih.gov/datasets/genome/.
2. Type the species of interest in the Genome search bar, e.g, Electrophorus electricus.
3. Look for the recent assembly that has RefSeq data as well.
4. Click on Action > Download
5. Check the boxes:
    - Annotation features (GTF)
    - Genomic coding sequences (FASTA)
    - Transcripts (FASTA) [*optional*]
6. Download the zip file.

I copied the 5 zip archives to below directory and unzipped them.

```bash
## pwd = Botany563/data/
mkdir -p ncbi_genome_downloads/zipped
cd ncbi_genome_downloads

for f in zipped/*.zip; do
    base=$(basename "$f" .zip)
    unzip -q "$f" -d "$base"
done
```

### Download meta-data

I downloaded all the supplementary materials from the Gallant et. al. (2014) and deposited them in:

```bash
## data/gallant_etal2014_metadata/
ls gallant_etal2014_metadata/
NIHMS884961-supplement-Supplementary_Information.pdf  NIHMS884961-supplement-Table_S1.xlsx  NIHMS884961-supplement-Table_S2.xls  paper.pdf
```

## PRE-PROCESSING

The first step is to split the metadata xlsx sheets into individual species-specific metadata files.

```bash
conda install -c conda-forge numpy pandas openpyxl -y
cd scripts/data_process
python split_metadata.py
```

### Selecting a representative subset of EO-associated genes

The original *E. electricus* metadata contains a few hundred genes across clusters 1 and 9, which is more than needed for an initial comparative CDS analysis. To create a manageable and biologically meaningful subset, I selected **25 genes from cluster 9** and **25 genes from cluster 1**, giving a total of **50 genes** for downstream analysis.

The selection was based on **electric-organ versus skeletal-muscle contrast**, rather than raw EO expression alone. This is important because a gene can have high EO expression but still not be EO-specific if it is also highly expressed in skeletal muscle. To capture this contrast, I computed the mean expression across the three EO tissues (`main EO`, `Sachs' EO`, and `Hunter's EO`) and compared it to `sk. muscle` using a log-ratio score.

- **Cluster 9 genes** were ranked by the highest **EO-to-muscle** contrast, since this cluster represents genes upregulated in electric organ.
- **Cluster 1 genes** were ranked by the strongest **muscle-to-EO** contrast, since this cluster represents genes downregulated in electric organ.
- To avoid repeated selection of the same gene, rows were **deduplicated by gene name** (`match_name`), keeping the best-scoring representative entry for each gene.
- The script then outputs the **top 25 genes from each cluster** and combines them into one final table for downstream sequence retrieval and alignment.

```bash
## cd Botany563/scripts/data-process
python get_top_EO_genes.py

```

I plan to use these genes of interest to get the corresponding sequences across the species.
