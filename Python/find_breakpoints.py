"""
This program browses .net files from pairwise whole genome alignments (http://hgdownload.cse.ucsc.edu/downloads.html) in order to find evolutionary breakpoints between the aligned species and stores the breakpoint coordinates in .bed files. Fills labelled 'top' are from the top level of the net and in combination with their lower class fills of types 'inv' and 'syn' (-/+ directionality compared to parent fill and from same chromosome) constitute a syntenic region. Lower level 'nonSyn' fills come from different chromosomes and are recognized as additional breakpoints. 
INPUT: Pairwise whole genome alignment in .net format from UCSC (Kent et al. 2003): human genome vs species x. Can be run from command line with the first param being the abbreviation of the first species, the 2nd param being that of the 2nd species .
OUTPUT: Store breakpoints in human genome (in comparison to genome of species x): Breakpoint coordinates are stored with their respective chromosome into BED files with the information in each line being chromosome, breakpoint start, breakpoint end (the last two being the same). 'top' level fills start points represent a breakpoint and also their end points. 'inv' and 'syn' fills are ignored but 'nonSyn' fills come from different chromosomes and represent another pair 
of breakpoints. The number of rearrangement blocks and their sizes are also stored in a list.
"""

import argparse

# size_threshold: only fills with sizes >threshold will be considered
# bp_list (returned): list of all bp, needed to calculate sizes of rearrangement blocks (synteny blocks)
def find_breakpoints(infile, outfile, size_threshold):
    
    print('INFO: Net alignment infile: ', infile)
    print('INFO: Writing breakpoints to: ', outfile)
    
    infile = open(infile, 'r')
    outfile = open(outfile, 'w')
    bp_list = []

    for line in infile:
        line = line.lstrip()
        if line.startswith('net'):
            line = line.split(" ")
            chrom = line[1]
            chr_len = int(line[2])
        elif line.startswith('fill'):
            line = line.split(" ")
            fill_type = line[line.index('type') + 1]
		
            if fill_type in ['top', 'nonSyn'] and int(line[12]) > size_threshold:
                
                # extract breakpoints
                fst_break_start = line[1]
                fst_break_end = str(int(line[1]) + 1)
                
                # test if is chrom start
                if int(fst_break_start) != 0 :
                    outfile.write(chrom + "\t" + fst_break_start + "\t" + fst_break_end + "\n")

                sec_break_start = str(int(line[1]) + int(line[2]))
                sec_break_end = str(int(line[1]) + int(line[2]) + 1)
                
                if sec_break_end != chr_len:
                    outfile.write(chrom + "\t" + sec_break_start + "\t" + sec_break_end + "\n")
                
                # for calculating syn region sizes
                bp_list.extend([(chrom, int(fst_break_start)), (chrom, int(sec_break_start))])

    outfile.close()
    infile.close()


def commandline():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Path to input file in  Net format. See https://genome.ucsc.edu/goldenpath/help/net.html")
    parser.add_argument("-t", "--size_threashold", type=int, required=False, default = 10000, help="Minimum size for fills.")
    parser.add_argument("-o", "--output_file", type=str, required=False, default="outfile.bp.bed", help="Filename of breakpoints in BED file.")
    return parser.parse_args()


if __name__ == '__main__':

    # read commandline arguments
    args = commandline()

    find_breakpoints(args.input_file, args.output_file, args.size_threashold)
    
