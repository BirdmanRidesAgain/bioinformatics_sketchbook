#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Import libs
from Bio import SeqIO
from argparse import ArgumentParser
from mimetypes import guess_type
from functools import partial
import gzip
import pandas as pd
#import itertools

exampleString="""
The following is an example of what the -t tsv is expected to look like.
The top line is a header, the remaining data is data.

species	genome	ctg_name	length	span_start	span_end
Chlamydotis macqueenii	fasta/houbaraONT_assembly-racon-haplotigs00000000.fsa	contig_4330	6656	595	643
"""


def open_fa(fa_file):
    '''
    Takes in a path to a fa file (compressed or uncompressed).
    Returns a list of all records in the file.
    '''
    record_lst = []

    encoding = guess_type(fa_file)[1]
    _open = partial(gzip.open, mode = 'rt') if encoding == 'gzip' else open
    with _open(fa_file) as f:
        genome_iterator = SeqIO.parse(f, 'fasta')
        for record in genome_iterator:
            record_lst.append(record)

    return record_lst

def splice_out(record, edit):
    """
    Function takes a seq. record object and a tuple of integers, where the first value is lower than the record..
    It removes specified regions of a sequence and returns the fixed record.
    """
    expected_length=(int(len(record.seq)) - (edit['span_end'] - edit['span_start']))

    first_half_seq=record.seq[0:edit['span_start'].iloc[0]]
    second_half_seq=record.seq[edit['span_end'].iloc[0]:edit['length'].iloc[0]]
    record.seq=first_half_seq+second_half_seq
    
    length_after_trim=len(record.seq)
    assert (int(expected_length.iloc[0])==length_after_trim), f"Trimming failed for {record.name}.\n\tExpected length: {int(expected_length)}\n\tActual length: {length_after_trim}"

    return(record)

def make_fa_name(species, suffix):
    '''
    Makes a name from the species and file suffix.
    '''
    return(f"{species.replace(' ','_')}_{suffix}")


def main():
    # 1. Use argparse to read in the path to the args
    parser = ArgumentParser(description="Takes in a fa and a Genbank-provided contamination file, and removes the specified contaminants.")
    parser.add_argument('-f','--fasta', help='Path to fasta file. You are expected to use the cleaned one GenBank gives you (eg,\'houbaraONT_assembly-racon-haplotigs00000000.fsa\').', type=str, required=True)
    parser.add_argument('-t', '--tsv', help='TSV of all the edits GenBank needs', type=str, required=True)
    parser.add_argument('-s', '--species', help='The species your\'e editing. Needed because GenBank alters the contig names.', type=str, required=True)
    args = parser.parse_args()

    # Read in inputs

    # NOTE: this code assumes that all species are unique. This is not guaranteed.
    df = pd.read_csv(args.tsv, sep="\t")
    edits=df[df['species'] == args.species]

    # Sanity check the the beginning of your sequence edit is before the end
    assert (df['span_start'] < df['span_end']).all(), f"Your span_start should come before your span_end."

    # sanity check for whether its finding the edits
    if len(edits) == 0:
        raise ValueError(f"No edits were found for {args.species}. Did you mean: {list(df['species'])}?")

    fa = open_fa(args.fasta)

    # Make the edits.
    # NOTE: this code assumes there is only one edit needed per contig, which is not guaranteed for all cases
    # Fix it later if needed
    counter=0

    for i in fa:
        if (i.name in list('lcl|'+edits['ctg_name'])):
            edit=edits[edits['ctg_name']==i.name[4:]]
            i=splice_out(i, edit)
            counter=counter+1
    
    # Sanity check to see if the number of trims matches the number of edits
    if (counter<len(edits)):
        raise ValueError(f"Fewer trims than expected were made. {counter} were made; {str(len(edits))} were expected.")
    else:
        print(f"{counter} trims were made; {str(len(edits))} were expected.")

    # now write the fasta
    with open(make_fa_name(args.species,'GenBank_edits.fa'), "w") as output_handle:
        SeqIO.write(fa, output_handle, 'fasta')
    print("done")
    return 0

# If this is being imported
if __name__=="__main__":
    main()
