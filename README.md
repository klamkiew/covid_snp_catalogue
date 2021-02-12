# Covid SNP Catalogue
[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/klamkiewicz?label=%40klamkiewicz&style=social)](https://twitter.com/klamkiewicz)
***

A small, rather hacky, script that takes a bunch of SARS-CoV-2 genomes and compares them to the RefSeq genome.

Dependencies:
```bash
conda create -n covid_snp_catalogue -c bioconda python=3.7 biopython mafft && conda activate covid_snp_catalogue
```

In theory, the `createDatabase.sh` script should do everything for you. 
The only input needed is a multiple fasta file with SARS-CoV-2 genomes of interest.

```bash
bin/createDatabase.sh <FASTA-FILE>
```

If you want to re-use the information gathered of an earlier run, simply add `knownSNPs.pkl` to the bash command:


```bash
bin/createDatabase.sh <FASTA-FILE> knownSNPs.pkl
```

WARNING: this currently breaks the .csv file produced by the script (see [this issue](https://github.com/klamkiew/covid_snp_catalogue/issues/1)).


***

#### Example Output

```
NEW ISOLATE
EPI_ISL_466918

Mutations:
RNA: c241t -- Protein: Mutation not part of CDS
RNA: t2321c -- Protein: ORF1ab: Synonymous mutation L686L
RNA: c3037t -- Protein: ORF1ab: Synonymous mutation Y925Y
RNA: c3092t -- Protein: ORF1ab: P943S
RNA: c14408t -- Protein: Mutation not part of CDS
RNA: g15327t -- Protein: Mutation not part of CDS
RNA: a23403g -- Protein: S: D614G
RNA: c27999t -- Protein: ORF8: P36S

Mutations are also observed in:
c241t: 239 out of 285 Isolates
t2321c: 0 out of 285 Isolates
c3037t: 228 out of 285 Isolates
c3092t: 12 out of 285 Isolates
c14408t: 230 out of 285 Isolates
g15327t: 12 out of 285 Isolates
a23403g: 243 out of 285 Isolates
c27999t: 12 out of 285 Isolates

```
