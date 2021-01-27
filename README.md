# Covid SNP Catalogue
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


***

####Example Output

```
NEW ISOLATE
EPI_ISL_466915

Mutations:
RNA: c241t -- Protein: Mutation not part of CDS
RNA: c3037t -- Protein: orf1a: Synonymous mutation Y925Y
RNA: c3092t -- Protein: orf1a: P943L
RNA: c14408t -- Protein: orf1ab: Synonymous mutation P314P
RNA: g15327t -- Protein: orf1ab: L621F
RNA: a23403g -- Protein: S: D614E
RNA: g25855t -- Protein: ORF3a: D155V
RNA: g26840a -- Protein: M: R107S
RNA: c27999t -- Protein: ORF8: P36L

Mutations are also observed in:
c241t: 239 out of 285 Isolates
c3037t: 228 out of 285 Isolates
c3092t: 12 out of 285 Isolates
c14408t: 230 out of 285 Isolates
g15327t: 12 out of 285 Isolates
a23403g: 243 out of 285 Isolates
g25855t: 7 out of 285 Isolates
g26840a: 11 out of 285 Isolates
c27999t: 12 out of 285 Isolates


```