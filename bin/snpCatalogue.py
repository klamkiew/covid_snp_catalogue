#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys
import pickle
from collections import defaultdict
from collections import Counter
from Bio import AlignIO

codon2aminoacid = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
  }

annotationFile = sys.argv[1]
alignmentPath = sys.argv[2]

def translateCDS(sequence):
  protein = ''
  for x in range(0,len(sequence),3):
    codon = sequence[x:x+3]
    if len(codon) != 3:
      exit(1)
    try:
      protein += codon2aminoacid[codon]
    except KeyError:
      protein += 'X'
  return(protein)

proteinRegions = {}
with open(annotationFile, 'r') as inputStream:
  for line in inputStream:
    if line.startswith("#") or not line:
      continue
    line = line.strip()
    currentLine = line.split()
    if currentLine[2] == "CDS":
      start = int(currentLine[3])-1
      stop = int(currentLine[4])
      geneName = line.split("gene=")[1].split(';')[0]
      proteinRegions[geneName] = (start,stop)


REFERENCE = "NC_045512.2_Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan-Hu-1,_complete_genome"

if len(sys.argv) == 4:
  knownSNPs = pickle.load(open(sys.argv[3], "rb"))
  pickeldSNPs = sys.argv[3]
else:
  knownSNPs = defaultdict(list)
  pickeldSNPs = 'knownSNPs.pkl'

alignment = AlignIO.read(alignmentPath, "fasta")
referenceAln = [str(x.seq) for x in alignment if x.id == REFERENCE][0]
referenceSequence = ''.join([x for x in referenceAln if x != '-']).upper()

proteins = {name : translateCDS(referenceSequence[x:y]) for name,(x,y) in proteinRegions.items()}
recordsInAln = []
newRecords = []

for record in alignment:
  recordsInAln.append(record.id)
  if record.id in knownSNPs:
    continue
  if record.id == REFERENCE:
    continue
  newRecords.append(record.id)
  for idx, nucleotide in enumerate(str(record.seq)):
    if nucleotide == 'n' or nucleotide == '-': continue
    if nucleotide != referenceAln[idx]:
      position = len(referenceAln[:idx]) - Counter(referenceAln[:idx])['-']
      if position <= 0: continue
      key = f"{referenceAln[idx]}{position+1}{nucleotide}"
      knownSNPs[key].append(record.id)
      knownSNPs[record.id].append(key)

for record in recordsInAln:
  if not knownSNPs[record]:
    knownSNPs[record] = []
pickle.dump(knownSNPs, open(pickeldSNPs ,'wb'))

if not newRecords:
  print("No new records found! Exiting...")
  exit(3)

def is_int(n):
  try:
    int(n)
  except ValueError:
    return False
  else:
    return True

mutation2isolates = {k:tuple(v) for k,v in knownSNPs.items() if is_int(k[1:-1])}
isolate2mutations = {k:tuple(v) for k,v in knownSNPs.items() if k in newRecords}

try:
  del mutation2isolates[REFERENCE]
except KeyError:
  pass

def getDDL():
  return defaultdict(str)
proteinMutation = defaultdict(getDDL)

for isolate, mutations in isolate2mutations.items():
  for mutation in mutations:
    pos = int(mutation[1:-1])
    wt = mutation[0]
    mut = mutation[-1]
    cds = [(x,y,z) for x,(y,z) in proteinRegions.items() if y <= pos <= z]
    if cds:
      cds = cds[0]

      mutatedSeq = translateCDS(referenceSequence[cds[1]:pos] + mut.upper() + referenceSequence[pos+1:cds[2]])
      refSeq = proteins[cds[0]]
      for idx, aa in enumerate(refSeq):
        if aa != mutatedSeq[idx] and mutatedSeq[idx] != "X":
          proteinMutation[isolate][mutation] = f"{cds[0]}: {aa}{idx+1}{mutatedSeq[idx]}"
      if proteinMutation[isolate][mutation] == '':
        ntInOrf = int((pos - cds[1])/3)
        if mutatedSeq[ntInOrf] != "X":
          proteinMutation[isolate][mutation] = f"{cds[0]}: Synonymous mutation {refSeq[ntInOrf]}{ntInOrf+1}{mutatedSeq[ntInOrf]}"
        else: proteinMutation[isolate][mutation] = f"{cds[0]}: Undetermined AminoAcid due to sequencing uncertainty"
    else:
      proteinMutation[isolate][mutation] = f"Mutation not part of CDS"


for isolate, mutations in isolate2mutations.items():
  if not mutations: continue
  sanePath = isolate.replace('/','-').replace('|','-')
  with open(f"SUMMARY_{sanePath}.txt", 'w') as outputStream:
    if isolate in newRecords:
      outputStream.write(f"NEW ISOLATE\n{isolate}\n\nMutations:\n")
      for mutation in mutations:
        pos = int(mutation[1:-1])
        outputStream.write(f"RNA: {mutation}")
        if proteinMutation[isolate][mutation]:
          outputStream.write(f" -- Protein: {proteinMutation[isolate][mutation]}\n")
        else:
          outputStream.write("\n")
      outputStream.write("\nMutations are also observed in:\n")
      for mutation in mutations:
        otherIsolates = [iso for iso in mutation2isolates[mutation] if iso != isolate]
        outputStream.write(f'{mutation}: {len(otherIsolates)} out of {len(isolate2mutations)-1} Isolates\n')

sortedMutation = sorted(list(mutation2isolates), key=lambda x:int(x[1:-1]))
with open("mutation_catalogue.csv", 'w') as outputStream:
  outputStream.write(";"+";".join(sortedMutation)+'\n')
  for record in alignment:
    isolate = record.id
    if isolate == REFERENCE: continue
    mutations = isolate2mutations[isolate]
    outputStream.write(isolate+";"+";".join(["X" if x in mutations else "-" for x in sortedMutation])+'\n')
