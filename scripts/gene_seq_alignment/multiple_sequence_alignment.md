# Prepare the gene sequence file and run MSA

## Prepare the input files to the MSA algorithms

The multiple sequence alignment (MSA) algorithms require inputs in the form of a FASTA file. In my case, the inputs would be per-gene sequence files which have the corresponding qequences for all species


### Get the per-species sequences for common marker genes

The script uses the ranked electric-eel marker lists from cluster 9 and cluster 1 as the reference gene set. It scans each species CDS FASTA, extracts gene names from the FASTA headers, and checks which of the reference markers are present in each species. Gene matching is done case-insensitively so that differences in capitalization across annotation files do not affect the comparison. The common gene set is then defined as the intersection of these markers across all selected species.

After finding the shared markers, the script writes them out in the same ranked order as the original up and down reference lists and keeps the corresponding up/down label for each gene. It also creates a filtered CDS FASTA for each species containing only sequences from these intersecting genes. From this shared set, I plan to use the top 25 upregulated and top 25 downregulated markers to generate the final input for the MSA algorithms.

```bash
# cd scripts/gene_seq_alignment
python generate_filtered_species_fa_files.py
```
SCRIPT OUTPUT:
```text
> Bb: 89 EO markers found
> Bg: 134 EO markers found
> Ee: 134 EO markers found
> Pk: 86 EO markers found
> Ip: 130 EO markers found

> Genes present in ALL species: 72
> Bb: wrote 185 CDS records -> ../../data/filtered_per_species_CDS_fa/Bb_intersect_genes_cds.fna
> Bg: wrote 176 CDS records -> ../../data/filtered_per_species_CDS_fa/Bg_intersect_genes_cds.fna
> Ee: wrote 269 CDS records -> ../../data/filtered_per_species_CDS_fa/Ee_intersect_genes_cds.fna
> Pk: wrote 223 CDS records -> ../../data/filtered_per_species_CDS_fa/Pk_intersect_genes_cds.fna
> Ip: wrote 192 CDS records -> ../../data/filtered_per_species_CDS_fa/Ip_intersect_genes_cds.fna
```

### Create the pre-gene sequence input files

Now that we have the information of common markers annotated across the species of interest, I plan to use top 25 up regulated and top 25 down regulated EO markers as my final genelist. As seen from the output above, it is clear that several genes have multiple CDS sequences which may be different isoforms of the same gene.

I deal with this by taking the longest sequence per gene for each species. The final output written will have 5 seqences per gene; 50 files and with headers containing the lcl, species, gene, transcript and protein information.

```bash
## pwd = scripts/gene_seq_alignment
python get_1-1cds_prep_input.py
```
SCRIPT OUTPUT:
```text
:
:
wrote ../../data/input_to_MSA_perGene_fa/atp1a2a.fna
wrote ../../data/input_to_MSA_perGene_fa/PCDH20.fna
wrote ../../data/input_to_MSA_perGene_fa/sh2d5.fna
:
:
```

## MSA Alignment - MAFFT and MACSE

For the project, I have 50 protein-coding nucleotide sequences (CDS) with 5 species each. Each file contains one representative CDS per species after taking the longest CDS per gene. So, this is a small gene-wise alignment query.

For this purpose, **MAFFT and MACSE** seem the most appropriate alignment tools. This choice was guided by both the biological nature of the dataset and the practical constraints of the analysis. Because these are coding sequences, a codon-aware method is especially relevant: MACSE was designed specifically for protein-coding nucleotide data and explicitly accounts for codon structure, frameshifts, and stop codons, making it well suited for evolutionary analyses of CDS.

In contrast, MAFFT is a widely used general-purpose multiple sequence alignment program, and its accuracy-oriented L-INS-i strategy is recommended for relatively small datasets such as mine. Using these two together allows me to compare a CDS-specialized alignment method against a well-established general MSA method

| Software | Description | Strengths | Weaknesses | Assumptions | User choices |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **MACSE** | A multiple sequence alignment method specifically designed for **protein-coding nucleotide sequences**. It aligns CDS using their amino acid translations while preserving codon structure and accounting for frameshifts and stop codons. | Well matched to **coding-sequence** data; preserves reading frame; can detect or tolerate frameshifts and stop codons; especially useful when downstream analyses depend on codon structure. | More specialized than general-purpose aligners; can be less straightforward to use and interpret than standard MSA tools; may be unnecessary if sequences are already very clean and perfectly in frame. | Assumes the input sequences are **homologous protein-coding sequences** and that translation-aware alignment is biologically meaningful for the dataset. | Choice of subprogram (here `alignSequences`), handling of output nucleotide and amino-acid alignments, and whether to use MACSE as the primary alignment or mainly as a CDS-aware comparison method. |
| **MAFFT** | A general-purpose multiple sequence alignment program for nucleotide or amino acid sequences. For this project, the **L-INS-i** strategy is appropriate because it is an accuracy-oriented method recommended for smaller datasets. | Widely used, well documented, and easy to run; strong overall performance; **L-INS-i** is accuracy-oriented and appropriate for small alignments such as 5 sequences per gene. | Not specifically codon-aware; may not preserve reading frame in the same explicit way as MACSE; like all automatic aligners, it can still produce questionable regions that require inspection. | Assumes the sequences are homologous and alignable under a standard nucleotide alignment framework; **L-INS-i** is most appropriate when local pairwise alignment information improves accuracy. | Choice of MAFFT strategy is important; here, `--localpair --maxiterate 1000` (L-INS-i) was chosen instead of a faster large-dataset mode because the dataset is small and accuracy is preferred over speed. |

### MAFFT ALIGNMENT

Useful links/resources:
1. https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html
2. https://mafft.cbrc.jp/alignment/software/tips0.html

```bash
## pwd = Botany563/scripts/gene_seq_alignment
mafft --version
```
```text
v7.526 (2024/Apr/26)
```
```bash
chmod +x run_mafft.sh
./run_mafft.sh
```

MAFFT v7.526 was run with the L-INS-i strategy `(--localpair --maxiterate 1000)` on each per-gene FASTA file. Log output showed successful convergence and completion for the loci examined, with MAFFT reporting the L-INS-i local-pair iterative-refinement strategy for these nucleotide CDS alignments.

MAFFT alignments completed successfully and produced reasonable results for several loci. However, a random inspection of few CDS alignments showed that some loci contained non-triplet gap patterns, *indicating that codon structure was not always preserved*. Therefore, MAFFT was retained as a comparison method, while MACSE was used as a coding-sequence-aware alignment method for the primary CDS analysis.

### MACSE ALIGNMENT

Useful links/resources:
1. https://www.agap-ge2pop.org/macse/
	- specifically: https://www.agap-ge2pop.org/macse/pipeline-documentation/
2. https://academic.oup.com/mbe/article/35/10/2582/5079334

```bash
## pwd = Botany563/scripts/gene_seq_alignment
macse -help
```
```text
This is MACSE V2.07 If you find MACSE useful, please cite:

MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons
Vincent Ranwez, Emmanuel J P Douzery, Cedric Cambon, Nathalie Chantret, Frederic Delsuc
Molecular Biology and Evolution, 2020, 35(10):2582-2584, https://doi.org/10.1093/molbev/msy159

MACSE: Multiple Alignment of Coding SEquences accounting for frameshifts and stop codons.
Vincent Ranwez, Sebastien Harispe, Frederic Delsuc, Emmanuel JP Douzery
PLoS One 2011, 6(9):e22594, https://doi.org/10.1371/journal.pone.0022594


usage: [-prog PROGRAM_NAME] [-debug] [-help]

PROGRAMS:
  'alignSequences'       aligns nucleotide (NT) coding sequences using their amino acid (AA) translations
  'alignTwoProfiles'     aligns two previously computed nucleotide alignments (also called profiles) without questioning them
  'enrichAlignment'      adds sequences to a pre-existing nucleotide alignment
  'exportAlignment'      allows to export a MACSE alignment and to compute some statistics, it can...
  'mergeTwoMasks'        indicates nucleotides kept after applying mask1 filtering then mask2 filtering (useful for traceability)....
  'multiPrograms'        sequentially executes multiple MACSE commands contained in a text file (one per line)....
  'refineAlignment'      improves the input nucleotide alignment
  'reportGapsAA2NT'      uses a amino acid alignment to align nucleotide sequences....
  'reportMaskAA2NT'      uses a filtered amino acid alignment to filter a nucleotide alignment....
  'splitAlignment'       splits one alignment, to extract a subset of sequences and/or sites....
  'translateNT2AA'       translates nucleotides into amino acids
  'trimAlignment'        trims the input alignment by removing gappy sites at the beginning/end of the alignment....
  'trimNonHomologousFragments'   identifies (and trims) sequence fragments that do not share homology with other sequences and remove those fragments....
  'trimSequences'        removes the 3' and 5' parts of the input sequence that are non homologous to an alignment....
```

```bash
chmod +x run_macse.sh
```

MACSE v2.07 was run on each per-gene FASTA file using the `alignSequences` program with the Standard genetic code. Log output showed successful completion across all loci, with MACSE producing both codon-aware nucleotide alignments and corresponding translated amino acid alignments for each gene. I chose MACSE because it is specifically designed for protein-coding nucleotide sequences and explicitly preserves codon structure.

A random inspection of several MACSE nucleotide alignments showed that gap patterns were consistent with coding-sequence alignment, with sampled loci preserving triplet structure more cleanly than the corresponding MAFFT nucleotide alignments. Although some loci remained gap-rich, the observed gap structure was codon-consistent in the examples examined.

Hence, I have decided to retain MACSE as the primary alignment method for downstream CDS-based analyses, while MAFFT was used as a general-purpose comparison alignment method.