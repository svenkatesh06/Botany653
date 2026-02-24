# PHYLOGENY PROJECT - 563 (Spring 2026)

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
conda create -p "<...path_to_env>/e_phylo" \
 python=3.11 \
 sra-tools \
 mafft muscle clustalw \
 r-base r-essentials \
 r-ape r-adegenet r-phangorn \
 bioconductor-ggtree

conda activate envs/e_phylo
```

### Data download

- I noticed that some of the species have single reads while the others have paired-end reads. 
- Additionally, 4 out of 5 sequences require SRA-toolkit(https://github.com/ncbi/sra-tools/wiki) to download the sequences as they are >5G bp. I have installed it in the above environment. 

```bash
# pwd = <path to git cloned repo-Botany653>/
mkdir scripts data results
cd data
mkdir raw_downloads logs tmp fastq qc
cd raw_downloads
vim srr_list.txt
```

I pasted the list of SRAs in the srr_list.txt
```bash
cat srr_list.txt 
SRR1299496
SRR1299090
SRR1299083
SRR1299088
SRR1299492
```
For getting data through SRA, I referred to the doc and ChatGPT. I have to use prefetch + fasterq dump.

**PREFETCH SRA**
```bash
# pwd = <path to git cloned repo-Botany653>/data/raw_downloads
while read -r srr; do \   
echo "==== prefetch $srr ====" | tee -a ../logs/prefetch.log;   prefetch --max-
size 200G "$srr" 2>&1 | tee -a ../logs/prefetch.log; \ 
done < srr_list.txt
```


Validate the pre-fetched SRRs and make sure they were downloaded completely.

```bash
# PROJ=<PATH TO CLONED REPO - Botany563>
while read -r srr; do
  echo "==== validate $srr ====" | tee -a "$PROJ/data/logs/validate.log"
  vdb-validate "$srr" 2>&1 | tee -a "$PROJ/data/logs/validate.log"
done < "$PROJ/data/raw_downloads/srr_list.txt"
```

**CONVERT SRA TO FASTQ**
```bash
# SCR = <PATH TO SCRATCH DIR>/<USER>/
mkdir -p "$SCR"/{fastq,tmp}

while read -r srr; do
  echo "==== fasterq-dump $srr ====" | tee -a "$PROJ/data/logs/fasterq.log"
  mkdir -p "$SCR/fastq/$srr"

  fasterq-dump "$srr" \
    --outdir "$SCR/fastq/$srr" \
    --split-files \
    --threads 8 \
    -t "$SCR/tmp" \
    --progress \
    2>&1 | tee -a "$PROJ/data/logs/fasterq.log"

  echo "==== compress $srr ====" | tee -a "$PROJ/data/logs/compress.log"
  pigz -p 8 "$SCR/fastq/$srr"/*.fastq 2>&1 | tee -a "$PROJ/data/logs/compress.log"

done < "$PROJ/data/raw_downloads/srr_list.txt"
```

Four of the SRA's have paired reads and 1 has single reads.
