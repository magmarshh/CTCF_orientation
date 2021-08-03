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


# Initialize the counts as global variables so multiple methods can access them
convergent = 0
tandem = 0
divergent = 0
single = 0
none = 0


def main():
    # parse command line arguments for loop, motif, and output files
    loopfile, motiffile, outputfile = initialize()

    # read in loop dataframe
    loop_df = read_loop(loopfile)

    # read in motif dataframe
    motif_df = read_motif(motiffile)

    # loop through each row in the loop dataframe and update the counts of each kind of motif
    loop_df.apply(lambda x: process_loop_row(motif_df, x), axis=1)
    ctcf_orientation_counts = [convergent, tandem, divergent, single, none]

    # plot and display the output of the orientation counts
    plot_output(ctcf_orientation_counts, outputfile)


# searches the ctcf motifs dataframe for overlaps with the inputted values
def find_strand(ctcf_df, loop_chrom, loop_row, strand_num):
    # subset the ctcf dataframe to only contain the correct chromosome for this row
    ctcf_chr = ctcf_df.loc[(ctcf_df['chrom'] == loop_chrom)]
    # if we're looking at anchor list 1...
    if strand_num == 1:
        # return the strand values for every location in the ctcf dataframe subset which contains EITHER:
        # a CTCF start site between the first start and end boundaries of the loop row OR
        # a CTCF end site betwen the first start and end boundaries of the loop row
        return ctcf_chr.loc[((ctcf_chr['start'] >= loop_row['start1']) & (ctcf_chr['start'] <= loop_row['end1'])) &
                            ((ctcf_chr['end'] >= loop_row['start1']) &
                             (ctcf_chr['end'] <= loop_row['end1'])), 'strand'].tolist()
    # same thing as above for anchor list 2,
    return ctcf_chr.loc[((ctcf_chr['start'] >= loop_row['start2']) & (ctcf_chr['start'] <= loop_row['end2'])) &
                        ((ctcf_chr['end'] >= loop_row['start2']) &
                         (ctcf_chr['end'] <= loop_row['end2'])), 'strand'].tolist()


# updates the counts of each type of CTCF orientation
def find_orientation(loop_anchor_1_list, loop_anchor_2_list):
    global none, single, convergent, tandem, divergent
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


# wrapper function for both find_orientation and find_strand
def process_loop_row(ctcf_df, loop_row):
    # Create empty lists to hold the strands for the CTCF motifs that match each anchor
    loop_anchor_1_list = []
    loop_anchor_2_list = []
    # Get the chromosomes of each anchor
    loop_anchor1_chr = loop_row['chrom1']
    loop_anchor2_chr = loop_row['chrom2']
    # find the strands for the CTCF motifs that match each anchor
    loop_anchor_1_list.extend(find_strand(ctcf_df, loop_anchor1_chr, loop_row, 1))
    loop_anchor_2_list.extend(find_strand(ctcf_df, loop_anchor2_chr, loop_row, 2))
    # update the counts of each type of CTCF orientation
    find_orientation(loop_anchor_1_list, loop_anchor_2_list)




# plots outputs of orientation counts
def plot_output(ctcf_orientation_counts, outputfile):
    labels = ['Convergent', 'Tandem', 'Divergent', 'Single', 'None']
    # Plot
    plt.clf()
    plt.pie(ctcf_orientation_counts, autopct='%1.1f%%',
            shadow=False, startangle=90)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.legend(labels)
    plt.savefig(outputfile)


# reads in motif file
def read_motif(motiffile):
    return pandas.read_csv(motiffile,
                           names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'pvalue', 'qvalue', 'seq'],
                           delimiter='\t')


# reads in loop file, removes header if it's present
def read_loop(loopfile):
    try:
        return pandas.read_csv(loopfile,
                           usecols=list(range(6)),  # use only the first 6 columns
                           names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'],
                           delimiter='\t',
                           dtype={'chrom1': str, 'start1': 'Int64', 'end1': 'Int64', 'chrom2': str,
                                  'start2': 'Int64',
                                  'end2': 'Int64'})
    except ValueError:
        return pandas.read_csv(loopfile,
                           usecols=list(range(6)),  # use only the first 6 columns
                           names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'],
                           delimiter='\t',
                           header=0,
                           dtype={'chrom1': str, 'start1': 'Int64', 'end1': 'Int64', 'chrom2': str,
                                  'start2': 'Int64',
                                  'end2': 'Int64'})

# parses command line inputs for loop, motif, and orientation files.
def initialize():
    loopfile = ''
    motiffile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hl:m:o:", ["lfile=", "mfile=", "ofile="])
    except getopt.GetoptError as err:
        print('ctcf_orientation.py -l <loopfile> -m <motiffile> -o <outputfile>')
        print(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('ctcf_orientation.py -l <loopfile> -m <motiffile> -o <outputfile>',
                  '\n **This script allows you to visualize the distribution of CTCF orientation within the loops of '
                  'the '
                  'given loop file. Outputs a piechart of CTCF orientation.**',
                  '\n **Requires: python 3.9, pandas 1.2.4 and matplotlib**',
                  '\n **Required Arguments:**', '\n \t -l <loopfile> : Pairedbed file of loops with the format '
                                                'chromosome anchor1_start anchor1_end chromosome anchor2_start '
                                                'anchor2_end',
                  '\n \t -m <motiffile> : Motif file in bed format, must be: chromosome start end name strength '
                  'orientation '
                  'pvalue qvalue sequence',
                  '\n \t -o <outputfile> : Desired name of output pie chart in matplotlib accepted picture format.')
            sys.exit()
        elif opt in ("-l", "--lfile"):
            loopfile = arg
        elif opt in ("-m", "--mfile"):
            motiffile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print('Loop file is ', loopfile)
    print('Motif file is ', motiffile)
    print('Output file is ', outputfile)
    return loopfile, motiffile, outputfile


if __name__ == "__main__":
    main()
