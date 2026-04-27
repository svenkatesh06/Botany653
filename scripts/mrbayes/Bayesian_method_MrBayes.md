# BAYESIAN INFERENCE WITH MrBayes

Useful Link:
1. https://github.com/PatrickKueck/FASconCAT-G/tree/master
2. https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf

## PREPARE INPUT TO MRBAYES : FASconCAT-G

For the Bayesian phylogenetic analysis, I prepared a separate MrBayes-compatible input file from the cleaned nucleotide alignments. MrBayes requires an aligned sequence matrix in **NEXUS** format.

**NOTE:**

I first attempted to concatenate the cleaned nucleotide alignments using [AMAS](https://github.com/marekborowiec/AMAS/blob/master/README.md), but the interleaved NEXUS output was not parsed correctly by MrBayes in this case, *so I switched to FASconCAT-G for preparing the final supermatrix*.
FASconCAT-G was used to concatenate the cleaned per-gene nucleotide alignments into a single supermatrix for Bayesian phylogenetic analysis.

```bash
# pwd = scripts/mrbayes
mkdir -p ../../results/mrbayes/fasconcat_results
cd ../../results/mrbayes/fasconcat_results

# download FASconCAT-G
wget https://raw.githubusercontent.com/PatrickKueck/FASconCAT-G/master/FASconCAT-G_v1.06.1.pl
chmod +x FASconCAT-G_v1.06.1.pl

# make .fas symlinks from the cleaned .fna files
for f in ../../../results/msa_macse/reordered_cleaned_macse_nt/*.fna; do
    base=$(basename "$f" .fna)
    ln -s "$f" "${base}.fas"
done
```


**CONCATENATING AND INTERLEAVING ALIGNMENTS**

The input files were the MACSE-cleaned nucleotide alignments from `results/msa_macse/reordered_cleaned_macse_nt/`. Because these alignments represent coding sequences and were already codon-aware from the MACSE step, they were suitable for downstream codon-position partitioning in MrBayes.

 Since FASconCAT-G recognizes standard FASTA-style file extensions more reliably, the `.fna` inputs were linked with `.fas` extensions in the working directory.

```bash
# pwd = results/mrbayes/fasconcat_results
perl FASconCAT-G_v1.06.1.pl -s -n -l > _fasconcat_test.log 2>&1

cp FcC_supermatrix.nex ../test_run_mrbayes.nex
```

The options used were:
- `-s`: starts the FASconCAT-G concatenation run.
- `-n`: produces NEXUS output for the concatenated supermatrix.
- `-l`: writes partition information for the concatenated alignment.

The **concatenated alignment contained 5 taxa and 107,499 nucleotide sites**. FASconCAT-G also produced gene-level partition coordinates, recording the start and end positions of each gene in the supermatrix.



**CODON-POSITION PARTITIONING:**

Because the NT alignments represent protein-coding sequences, I did not plan to use each gene as a fully independent partition for the main Bayesian analysis. For a dataset with only 5 taxa, assigning a separate substitution model to every gene could introduce too many parameters relative to the amount of information in the data. Instead, I prepared the mrbayes input using codon-position partitioning, which is more biologically appropriate for coding DNA. The mrbayes manual says that:
> *for protein-coding NT data, partitioning by codon position is often preferable to partitioning simply by gene, because third codon positions often evolve faster than first and second positions.*


## RUN MRBAYES

The FASconCATG-generated NEXUS file contains the concatenated alignment matrix, **and I appended a begin mrbayes; ... end; block to the same file to define the codon-position partitions, substitution model, and MCMC settings before running the file directly in mrbayes.**


```text
begin mrbayes;

    set autoclose=yes nowarnings=yes;

    charset codon1 = 1-.\3;
    charset codon2 = 2-.\3;
    charset codon3 = 3-.\3;

    partition codonpos = 3: codon1, codon2, codon3;
    set partition = codonpos;

    lset applyto=(all) nst=mixed rates=gamma;

	unlink statefreq=(all) shape=(all);
	prset applyto=(all) ratepr=variable;

    mcmcp ngen=1000000 samplefreq=1000 printfreq=1000 diagnfreq=5000 nruns=2 nchains=4 savebrlens=yes relburnin=yes burninfrac=0.25;

    mcmc;
    sump relburnin=yes burninfrac=0.25;
    sumt relburnin=yes burninfrac=0.25;

end;
quit;
```

**MrBayes block: command explanation**

- `begin mrbayes;`
  - Starts the MrBayes command block inside the NEXUS file.
  - This block tells MrBayes how to analyze the alignment after reading the sequence matrix.

- `charset codon1 = 1-.\3;`
  - Defines the first codon-position character set.
  - Selects alignment columns 1, 4, 7, 10, etc.
  - Used because the input nucleotide alignments are coding-sequence alignments.

- `charset codon2 = 2-.\3;`
  - Defines the second codon-position character set.
  - Selects alignment columns 2, 5, 8, 11, etc.
  - Allows second codon positions to be modeled separately.

- `charset codon3 = 3-.\3;`
  - Defines the third codon-position character set.
  - Selects alignment columns 3, 6, 9, 12, etc.
  - Third codon positions often evolve faster, so they should not be forced to share the same rate behavior as first and second positions.

- `partition codonpos = 3: codon1, codon2, codon3;`
  - Creates a partitioning scheme named `codonpos`.
  - The alignment is divided into three partitions corresponding to the three codon positions.
  - This is more appropriate for protein-coding nucleotide data than treating the entire concatenated alignment as one homogeneous block.

- `set partition = codonpos;`
  - Activates the `codonpos` partition scheme.
  - Without this command, the defined character sets would not necessarily be used as active model partitions.

- `lset applyto=(all) nst=mixed rates=gamma;`
  - Sets the substitution model for all active partitions.
  - `applyto=(all)` applies the model settings to all codon-position partitions.
  - `nst=mixed` allows MrBayes to sample across reversible nucleotide substitution models instead of fixing a single model such as GTR.
  - `rates=gamma` allows sites to evolve at different rates using a gamma-distributed rate model.
  - I used `gamma` rather than `invgamma` to avoid adding an extra invariant-sites parameter, which may be difficult to estimate reliably with only five taxa.

- `unlink statefreq=(all) shape=(all);`
  - Allows each codon-position partition to have its own parameters.
  - `statefreq=(all)` gives each partition its own nucleotide base frequencies.
  - `shape=(all)` gives each partition its own gamma shape parameter for among-site rate variation.
  - This is useful because codon positions can differ in base composition, and rate heterogeneity.

- `prset applyto=(all) ratepr=variable;`
  - Allows the overall evolutionary rate to vary among partitions.
  - This is important because third codon positions usually evolve faster than first and second codon positions.

- `mcmcp ngen=1000000 samplefreq=1000 printfreq=1000 diagnfreq=5000 nruns=2 nchains=4 savebrlens=yes;`
  - Sets the MCMC run parameters before starting the analysis.
  - `ngen=1000000` runs the analysis for 1,000,000 generations.
  - `samplefreq=1000` saves one sample every 1,000 generations.
  - `printfreq=1000` prints progress to the screen every 1,000 generations.
  - `diagnfreq=5000` calculates convergence diagnostics every 5,000 generations.
  - `nruns=2` runs two independent MCMC analyses, which helps assess convergence.
  - `nchains=4` runs one cold chain and three heated chains per run, improving exploration of tree space.
  - `savebrlens=yes` saves branch lengths in the sampled trees.

- `mcmc;`
  - Starts the Bayesian MCMC analysis using the settings defined above.

- `sump;`
  - Summarizes the sampled model parameters after the MCMC run.
  - Used to inspect estimates such as substitution parameters, base frequencies, rate variation, ESS, and PSRF.

- `sumt;`
  - Summarizes the sampled trees after burn-in.
  - Produces the Bayesian consensus tree with clade posterior probabilities and branch-length summaries.

- `end;`
  - Ends the MrBayes command block.

I partitioned the alignment by codon position **but avoided fully unlinking all substitution parameters across partitions in the initial analysis to reduce overparameterization**. Base frequencies and gamma shape parameters were allowed to differ among codon positions, and overall partition rates were allowed to vary, while the substitution-rate matrix was shared across partitions. This should provide a compromise between biological flexiblity and model stability for a 5-taxon dataset.

### TEST RUN OF MRBAYES

I first did a test run of MrBayes using the above settings on the concatenated nucleotide alignment to make sure everything was working correctly before running the full analysis. the test run was done for 10,000 generations, which is not long enough for convergence but allows me to check that the MCMC is running and producing output files without errors.

```bash
# pwd = scripts/mrbayes
bash wrapper_test_run_mrbayes.sh
```
The test run completed successfully in less than one minute, indicating that the input file and MrBayes block were formatted correctly and that the full run should be computationally manageable.

## MRBAYES FULL RUN

After confirming that the test run worked, I ran the full analysis for 1,000,000 generations using the same settings.

```bash
## pwd = scripts/mrbayes
cp ../../results/mrbayes/test_run_mrbayes.nex ../../results/mrbayes/full_run_mrbayes.nex
## change the mcmcp ngen=10000 to ngen=1000000 in full_run_mrbayes.nex before running
bash wrapper_run_full_mrbayes.sh
```

### MrBayes run results and convergence diagnostics

I used chatGPT to help interpret the MrBayes output and convergence diagnostics and make sense of all the parameters, which are summarized below:

- The full MrBayes analysis completed.
  - The run used 2 independent MCMC runs with four chains each.
  - The analysis was run for 1,000,000 generations.
  - Trees and parameters were sampled every 1,000 generations.
  - The full run completed in ~32 minutes.

- A relative burn-in of 25% was used.
  - This means the first 25% of samples were discarded before summarizing parameters and trees.

- The run showed good convergence between independent chains.
  - The average standard deviation of split frequencies was `0.000000`.
  - The maximum standard deviation of split frequencies was also `0.000000`.
  - This indicates that the two independent runs sampled the same posterior tree distribution.
  - The average PSRF was `1.000`.
  - The maximum PSRF was `1.001`.
  - **PSRF values close to 1.0 indicate that independent runs converged to similar posterior distributions.**

- The Bayesian consensus tree was summarized using `sumt`.
  - The consensus tree represents the posterior distribution of sampled trees after burn-in.
  - Clade support values in the MrBayes tree are posterior probabilities.
  - Posterior probabilities are not the same as bootstrap values from IQ-TREE.
  - In this run, the major clades had posterior probabilities of 1.00, indicating very strong posterior support.

- The recovered Bayesian topology was strongly supported.
  - The consensus tree supported a close relationship between `Bg` and `Ee`.
  - It also supported a close relationship between `Pk` and `Bb`.
  - With `Ip` treated as the outgroup, the Bayesian tree can be interpreted as supporting the topology `Ip | ((Pk, Bb), (Bg, Ee))`.
  - This topology can be compared directly with the IQ-TREE ML topology.


## ANALYZE MRBAYES OUTPUT USING TRACER

I downloaded the Tracer software for Windows OS:

https://github.com/beast-dev/tracer/releases/tag/v1.7.2 > Tracer.v1.7.2.zip

- Tracer was used to inspect the parameter trace files.
  - The `.p` files from both independent runs were loaded into Tracer.
  - The likelihood traces were stationary after burn-in.
  - The two independent runs showed overlapping trace behavior.
  - Most ESS values were above 200.
  - A few base-frequency parameters had lower ESS values, but they remained above 100.
  - Because the lowest ESS values were still above the common warning threshold and PSRF values were close to 1.0, the 1,000,000-generation run was considered sufficient for this analysis.

> *NOTE: In Tracer, the burn-in for each `.p` file was set manually to 250,000 generations to match the MrBayes run.*

- `LnL` represents the log likelihood.
  - It measures how well the sampled tree, branch lengths, and substitution model explain the observed alignment.
  - The absolute value of `LnL` is not interpreted biologically.
  - It is mainly used to check whether the MCMC trace has reached a stable region.

- `LnPr` represents the log prior probability.
  - It is the prior probability of the sampled tree, branch lengths, and model parameters.
  - Like `LnL`, it is mainly used as a convergence diagnostic rather than as a biological result.

- `TL{all}` represents total tree length.
  - It is the total amount of inferred evolutionary change across the tree.
  - The trace for tree length mixed well across the two independent runs.

- The nucleotide substitution-rate parameters were shared across codon-position partitions.
  - These are the `r(A<->C)`, `r(A<->G)`, `r(A<->T)`, `r(C<->G)`, `r(C<->T)`, and `r(G<->T)` parameters.
  - Because the substitution-rate matrix was not unlinked across partitions, these parameters are labeled `{all}`.
  - These parameters describe the relative rates of different nucleotide substitutions.
  - Their main use here was to check model behavior and convergence.

- Because `nst=mixed` was used, MrBayes sampled across reversible nucleotide substitution models.
  - This avoids forcing one fixed substitution model, such as GTR, before the analysis.
  - The sampled substitution-model parameter summaries indicate which reversible models had higher posterior support.
  - This was used as model averaging over nucleotide substitution models.

- Base frequencies were estimated separately for each codon-position partition.
  - Parameters such as `pi(A){1}`, `pi(C){1}`, `pi(G){1}`, and `pi(T){1}` correspond to codon position 1.
  - The `{2}` and `{3}` versions correspond to codon positions 2 and 3.
  - This was appropriate because coding sequences can have different nucleotide compositions across codon positions.

- Gamma shape parameters were also estimated separately for each codon-position partition.
  - These are `alpha{1}`, `alpha{2}`, and `alpha{3}`.
  - The gamma shape parameter describes among-site rate variation within a partition.
  - **Lower alpha values indicate stronger rate heterogeneity among sites.**
  - **Codon positions 1 and 2 show stronger among-site rate heterogeneity, while codon position 3 has more even rates within that partition.**

- Partition rate multipliers were estimated for each codon position.
  - These are `m{1}`, `m{2}`, and `m{3}`.
  - **The estimated rates showed that codon position 3 evolved fastest, codon position 1 was intermediate, and codon position 2 evolved slowest.**
  - **This matches the expected pattern for protein-coding nucleotide sequences, where third codon positions are often less constrained.**