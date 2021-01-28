#!/usr/bin/env bash
# Author: Kevin Lamkiewicz
# E-Mail: kevin.lamkiewicz@uni-jena.de

# The following shell script needs a .fasta file
# as input.
# It will then create an MSA in gapped-fasta format.
# Further, a Python3 script is invoked that'll scan all sequences
# for SNPs with regard to the reference Isolate from Wuhan.
# A summary is reported for each isolate, and an ugly 
# table for an overview is created.


# definitions. Change if needed
WD="$(dirname $0)"
SEQ="$1"
DIR="$(dirname $(realpath $SEQ))"
CPU=8
REF="$WD"/../reference/ncbi_reference_nc045512.fasta

if [ $# -gt 1 ]; then
	SNPS="$2"
	if [ ! -f "$SNPS" ];then
		echo "Invalid path name. Please check your SNP catalogue file."
		echo "$SNPS"
		exit 2
	fi
fi



SEQ=$(basename $SEQ)
BN=${SEQ%.*}

# if there is no input file to begin with
# we'll exit as well.
if [ ! -s "$DIR"/"$SEQ" ]; then
	echo "Invalid path name. Please check your input."
	echo "$DIR"/"$SEQ"
	exit 2
fi

# we need to include the NCBI
# reference strain here as well for downstream SNP analyses
INPUT="$DIR"/"$BN"_ncbiRef.fasta
cat "$REF" "$DIR"/"$SEQ" > "$INPUT"
sed -i '/>/ s/ /_/g' "$INPUT"

#cd "$DIR"
# creating multiple sequence alignment and phylogentic tree
# based on whole genome sequence on nucleotide level

mafft --reorder --thread "$CPU" "$INPUT" > "$DIR"/"$BN"_ncbiRef_mafft_fasta.aln

if [ ! -d "$DIR/snp_summaries" ]; then
  mkdir -p "$DIR/snp_summaries"
fi

# SNP analyses and report the overall SNP overview
if [ "$SNPS" ]; then
	python3 "$WD"/snpCatalogue.py "$WD"/../reference/sars-cov-2_nc045512.gff3 "$DIR"/"$BN"_ncbiRef_mafft_fasta.aln "$DIR/snp_summaries" "$SNPS"
else
	python3 "$WD"/snpCatalogue.py "$WD"/../reference/sars-cov-2_nc045512.gff3 "$DIR"/"$BN"_ncbiRef_mafft_fasta.aln "$DIR/snp_summaries"
fi

if [ $? -eq 3 ]; then
  exit 3
fi
