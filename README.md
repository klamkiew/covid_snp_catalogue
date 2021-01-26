# covid_snp_catalogue
[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/klamkiewicz?label=%40klamkiewicz&style=social)](https://twitter.com/klamkiewicz)
***

A small, rather hacky, script that takes a bunch of SARS-CoV-2 genomes and compares them to the RefSeq genome.

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