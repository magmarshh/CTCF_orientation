#!/usr/local/bin/python3.7

"""
Executable python script that takes in a paired bed file (loops in the format of:
chromosome1 anchor1_start   anchor1_end chromosome2 anchor2_start   anchor2_end)
and CTCF motif file and outputs the number of loops in the paired bed file that have
convergent, tandem, divergent, single, and none CTCF motifs.

"""

import getopt
import pandas
import sys
import matplotlib.pyplot as plt


def original_compute_ctcf(loop_df, ctcf_df):
    convergent = 0
    tandem = 0
    divergent = 0
    single = 0
    none = 0
    # Iterate through each loop (row) in the loop dataframe
    for loop_index, loop_row in loop_df.iterrows():
        # Create empty lists to hold the strands for the CTCF motifs that match each anchor
        loop_anchor_1_list = []
        loop_anchor_2_list = []
        # Get the chromosomes of each anchor
        loop_anchor1_chr = loop_row['chrom1']
        # print('Printing grubert anchor1', anchor1_chr)
        loop_anchor2_chr = loop_row['chrom2']
        # print('Printing grubert anchor2', anchor2_chr)
        # Iterate through each row in the CTCF motif dataframe
        for ctcf_index, ctcf_row in ctcf_df.iterrows():
            # If the chromosome of the CTCF motif is the same as the first anchor's chromosome
            if ctcf_row['chrom'] == loop_anchor1_chr:
                # If the CTCF motif intersects with the first anchor, add the strand to the anchor's list
                if ctcf_row['start'] in range(loop_row['start1'], loop_row['end1']) or ctcf_row['end'] in \
                        range(loop_row['start1'], loop_row['end1']):
                    loop_anchor_1_list.append(ctcf_row['strand'])
            if ctcf_row['chrom'] == loop_anchor2_chr:
                if ctcf_row['start'] in range(loop_row['start2'], loop_row['end2']) or ctcf_row['end'] in \
                        range(loop_row['start2'], loop_row['end2']):
                    loop_anchor_2_list.append(ctcf_row['strand'])
        # Now we have the CTCF orientations (strands) of each anchor in the loop, check and see if they are
        # First convergent, then tandem, then divergent, then single, then none
        # First check if either list is empty, automatically none
        # If the first anchor has no CTCF motifs (if empty list):
        # print('Printing anchor 1 list', loop_anchor_1_list)
        # print('Printing anchor 2 list', loop_anchor_2_list)
        if not loop_anchor_1_list:
            # Check if second anchor is empty also, if so then the loop is categorized as none
            if not loop_anchor_2_list:
                none += 1
            # If the second anchor is not empty, then the loop has a single CTCF motif
            if loop_anchor_2_list:
                single += 1
        # if the second anchor list is empty, already counted if they are both empty
        if not loop_anchor_2_list:
            # Check if the first anchor list is not empty, don't need to check if empty because already did above
            if loop_anchor_1_list:
                single += 1
        # Else, if they are both not empty, check if they are convergent, tandem, or divergent (in order)
        else:
            # Check if they are convergent, first will be + and - will be second
            if '+' in loop_anchor_1_list and '-' in loop_anchor_2_list:
                convergent += 1
            # Check if they are tandem
            elif '+' in loop_anchor_1_list and '+' in loop_anchor_2_list:
                tandem += 1
            elif '-' in loop_anchor_1_list and '-' in loop_anchor_2_list:
                tandem += 1
            # Check if they are divergent
            elif '-' in loop_anchor_1_list and '+' in loop_anchor_2_list:
                divergent += 1

    orientation_counts = [convergent, tandem, divergent, single, none]
    return orientation_counts


def main(argv):

    loopfile = ''
    motiffile = ''
    output_diff_file = ''
    output_breakdown_file = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hl:m:o:", ["lfile=", "mfile=", "ofile="])
    except getopt.GetoptError as err:
        print('ctcf_orientation.py -l <loopfile> -m <motiffile> -o <output_diff_file>')
        print(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('ctcf_orientation.py -l <loopfile> -m <motiffile> -o <output_diff_file>',
                  '\n **This script allows you to visualize the distribution of CTCF orientation within the loops of the '
                  'given loop file. Outputs a piechart of CTCF orientation.**',
                  '\n **Requires: python 3.9, pandas 1.2.4 and matplotlib**',
                  '\n **Required Arguments:**', '\n \t -l <loopfile> : Pairedbed file of loops with the format '
                                            'chromosome anchor1_start anchor1_end chromosome anchor2_start anchor2_end',
                  '\n \t -m <motiffile> : Motif file in bed format, must be: chromosome start end name strength '
                  'orientation '
                  'pvalue qvalue sequence',
                  '\n \t -o <output_diff_file> : Desired name of output pie chart in matplotlib accepted picture format.')
            sys.exit()
        elif opt in ("-l", "--lfile"):
            loopfile = arg
        elif opt in ("-m", "--mfile"):
            motiffile = arg
        elif opt in ("-o", "--ofile"):
            output_diff_file = arg
            temp = output_diff_file.split(".")
            output_breakdown_file = temp[0] + "_breakdown." + temp[1]
    print('Loop file is ', loopfile)
    print('Motif file is ', motiffile)
    print('Output file is ', output_diff_file)
    loop_df = pandas.read_csv(loopfile, names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'],
                                           delimiter='\t')
    motif_df = pandas.read_csv(motiffile,
                              names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'pvalue', 'qvalue', 'seq'],
                              delimiter='\t')
    ctcf_orientation_counts = original_compute_ctcf(loop_df,motif_df)
    ctcf_diff = [sum(ctcf_orientation_counts[:4]), ctcf_orientation_counts[4]]
    diff_labels = ['Differential', 'Non-differential']
    labels = ['Convergent', 'Tandem', 'Divergent', 'Single']

    ### Plot
    plt.clf()
    plt.pie(ctcf_orientation_counts[:4], autopct='%1.1f%%',
            shadow=True, startangle=90)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.legend(labels)
    plt.savefig(output_diff_file)
    plt.clf()

    plt.pie(ctcf_diff, autopct='%1.1f%%',
            shadow=True, startangle=90)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.legend(diff_labels)
    plt.savefig(output_breakdown_file)
    plt.clf()




if __name__ == "__main__":
    main(sys.argv[1:])

