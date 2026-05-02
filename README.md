# Phylogenetic Analysis of Electric Organ–Associated Genes in Electric Fishes

This repository contains my final project for **Botany 563: Phylogenetic Analysis of Molecular Data [Spring 2026]**.

The project uses electric-organ-associated coding sequences to compare phylogenetic signal across electric and non-electric fishes. The biological marker set was motivated by **Gallant et al. 2014**, which identified electric-organ-associated gene clusters, and the species/context choices were also guided by **Wang and Yang 2021**, which provides an updated comparative phylogenetic framework for electric fishes.

## Overview:

Electric organs have evolved independently across fishes and are useful for studying molecular convergence. In this project, I asked whether selected electric-organ-associated coding sequences recover the expected relationships among sampled electric fish lineages.

The expected broad topology is:

```text
Ip | ((Bb, Pk), (Bg, Ee))
```

where:

- `Bb` + `Pk` are African weakly electric fishes / Mormyroidea
- `Bg` + `Ee` are South American electric fishes / Gymnotiformes
- `Ip` is the non-electric outgroup

### Species included

| Code | Species | Group / role |
|---|---|---|
| `Ee` | *Electrophorus electricus* | Gymnotiformes; electric eel |
| `Bg` | *Brachyhypopomus gauderio* | Gymnotiformes |
| `Bb` | *Brienomyrus brachyistius* | Mormyroidea |
| `Pk` | *Paramormyrops kingsleyae* | Mormyroidea |
| `Ip` | *Ictalurus punctatus* | non-electric outgroup |


## Workflow summary

```text
Gallant et al. EO-associated marker genes
        |
        v
Select EO-associated markers
        |
        v
Match marker genes to current CDS annotations
        |
        v
Find genes present across all 5 species
        |
        v
Choose longest CDS per gene per species
        |
        v
Create per-gene FASTA files
        |
        v
Multiple sequence alignment
    - MAFFT
    - MACSE
        |
        v
Tree inference
    - Distance methods
    - Maximum parsimony
    - IQ-TREE maximum likelihood
    - MrBayes Bayesian inference
    - ASTRAL-IV coalescent species tree
        |
        v
Topology comparison
```

## Repository guide

Each major folder inside `scripts/` contains a corresponding `.md` file that describes the reproducible workflow for that step. These markdown files are the best place to start if you want to follow or rerun the analysis.

| Step | Folder to read | What it contains | Main output folder |
|---|---|---|---|
| 1. Data curation and marker selection | `scripts/data_process/` | EO marker filtering, CDS extraction, gene matching across species, longest-CDS selection | `data/filtered_per_species_CDS_fa/`, `data/input_to_MSA_perGene_fa/` |
| 2. Multiple sequence alignment | `scripts/gene_seq_alignment/` | MAFFT and MACSE alignment workflow; MACSE cleaning and reordering | `results/msa_mafft/`, `results/msa_macse/` |
| 3. Distance and parsimony analyses | `scripts/distance_parsimony/` | BioNJ distance trees, parsimony trees, consensus/tree comparison steps | `results/distance_parsimony/` |
| 4. Maximum-likelihood analysis | `scripts/ml_iqtree/` | IQ-TREE nucleotide and amino-acid analyses, partitioned ML analysis, tree visualization | `results/ml_iqtree/` |
| 5. Bayesian inference | `scripts/mrbayes/` | Concatenation, codon-position partitioning, MrBayes runs, Bayesian tree visualization | `results/mrbayes/` |
| 6. Coalescent species tree | `scripts/coalescent_astral/` | Per-gene IQ-TREE trees and ASTRAL-IV species-tree inference | `results/coalescent_astral/` |

