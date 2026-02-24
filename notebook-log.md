

# PHYLOGENY PROJECT

## DATA
The link to the paper chosen for this project is here:[Genomic basis for the convergent evolution of electric organs](https://pmc.ncbi.nlm.nih.gov/articles/PMC5541775/)

The goal of this paper is to uncover the genetic and transcriptional basis of convergent evolution in electric fishes—specifically, how a myogenic electric organ (a muscle-derived organ) has evolved repeatedly across distantly related lineages. To do this, the authors assembled the genome of the electric eel (*Electrophorus electricus*) and sequenced electric organ and skeletal muscle transcriptomes from multiple electric-fish species representing independent evolutionary origins of electrogenesis. 
I will focus on the **electric organ transcriptome sequences** across the sampled species and use them for comparative analyses aimed at understanding evolutionary relationships and shared molecular signatures associated with electric organ function - communication, navigation/electrolocation, and in some cases predation/defense.

Bioproject ID: PRJNA248545

SRAs of the electric organ sequencing:
1. https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR1299496
2. https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR1299090
3. https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR1299083
4. https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR1299088
5. https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR1299492


## ENVIRONMENT SETUP AND DATA DOWNLOAD

### Install the necessary packages for project in conda environment.

I initially tried installing conda environment in my windows system, but some of these packages require linux OS. 
So I tried on WSL Ubuntu, but that seemed to be pretty slow. 
Finally, I setup the environment on WID server. 

```bash
 conda create -p "/mnt/dv/wid/projects7/Roy-singlecell2/multispecies_singlecell/swathi_work/envs/e_phylo \
  python=3.11 \
  sra-tools \
  mafft muscle clustalw \
  r-base r-essentials \
  r-ape r-adegenet r-phangorn \
  bioconductor-ggtree
```

### Data download

- I noticed that some of the species have single reads while the others have paired-end reads. 
- Additionally, 3 of the sequences require SRA-toolkit(https://github.com/ncbi/sra-tools/wiki) to download the sequences as they are >5G bp. I have installed it in the above environment. 
