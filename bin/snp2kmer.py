#!/usr/bin/env python3

import os
import argparse
import textwrap
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq

PROGRAM = 'snp2kmer'

## Required functions

# Mutate sequence - OK
def mutateSequence(sequence, position, value):
    if position <= len(sequence):
        sequence = sequence[:position-1] + str(value) + sequence[position:]
    return sequence

# Get two possible sequences given a set of cordinates and two or more alleles
# - OK
def extractDnaSequences(ref_name, position, alleles):
    start = str(position - kmer_len + 1)
    end   = str(position + kmer_len)
    bedcoords = ' '.join([ref_name, start, end])
    
    reference_file = os.path.join(path_to_reference, "{}.fa".format(ref_name))

    bed   = pybedtools.BedTool(bedcoords, from_string=True)
    fasta = bed.sequence(fi=reference_file)
    fasta = SeqIO.read(fasta.seqfn, "fasta")
    base_sequence = str(fasta.seq).upper()

    variants = []
    for i in range(0, len(alleles)):
        var = mutateSequence(base_sequence, kmer_len, alleles[i])
        variants.append(var)
    variants.append(base_sequence)

    return tuple(variants)

# Determine which is the reference allele - OK
def findRefAllele(a1, a2, ref):
    if a1 == ref and a2 != ref:
        ref_allele = a1
        alt_allele = a2
    elif a1 != ref and a2 == ref:
        ref_allele = a2
        alt_allele = a1
    
    return (ref_allele, alt_allele)

# Convert a sequence to all its possible k-mers - OK
def getKmers(sequence):
    my_kmers = []
    for i in range(0, kmer_len):
        pos_end   = i + kmer_len
        substring = sequence[i:pos_end]
        my_kmers.append(substring)
    
    return my_kmers



def parse_arguments():
    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=textwrap.dedent("""

        Given a table containing SNPs with allelic imbalance and a table
        with k-mer f-VICE scores, calculate what is the f-VICE scores for
        both alleles regarding all the k-mers that span that coordinate.

        \0""")
    )

    parser.add_argument('-f', '--input_table', type=str, 
                        help='The table with all the allelic imbalance SNPS')
    parser.add_argument('-l', '--kmer_len', type=int, default=6, 
                        help='How big are your k-mers?')
    parser.add_argument('--reference_genome', type=str, 
                        default='/lab/data/reference/human/hg19', 
                        help='Path to your reference genome.')

    return parser.parse_args()


args = parse_arguments()


snp_file          = args.input_table
path_to_reference = args.reference_genome
kmer_len          = args.kmer_len


## Code 

# Read and process SNP file
with open(snp_file, 'r') as data_stream:

    # Update header
    header = data_stream.readline().strip().split()
    new_header = header + ["ref_allele", "alt_allele", "ref_kmer", "alt_kmer"]
    print(*new_header, sep="\t")

    # Process data
    for line in data_stream:
        line = line.strip().split("\t")
        
        chrom    = line[0]
        position = int(line[1])
        a1       = line[3]
        a2       = line[4]
        ref      = line[5]

        alleles = findRefAllele(a1, a2, ref)
        sequences = extractDnaSequences(chrom, position, alleles)
        ref_kmer_list = getKmers(sequences[0])
        alt_kmer_list = getKmers(sequences[1])

        for i in range(0, kmer_len):
            ref_kmer = ref_kmer_list[i]
            alt_kmer = alt_kmer_list[i]

            new_data = list(alleles) + [ref_kmer, alt_kmer]
            out_stream = line + new_data
            print(*out_stream, sep = "\t")

