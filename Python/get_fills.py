'''Creates BED file with all fill coordinates from a net file.
'''

import sys
import argparse

def commandline():
    parser = argparse.ArgumentParser(description=__doc__, 
    		formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", type=str, required=True,
                        help="Path to input file in  Net format. \
                        See https://genome.ucsc.edu/goldenpath/help/net.html")
    parser.add_argument("-o", "--output_file", type=str, 
    		required=False, default="outfile.bp.bed",
    		help="Filename of fills in BED file.")

    return parser.parse_args()

# size_threshold: only fills with sizes >threshold will be considered
# bp_list (returned): list of all bp, needed to calculate sizes of rearrangement 
# blocks (synteny blocks)
def checkFillSize(infile, outfile):
	
    print('INFO: Net alignment infile: ', infile)
    print('INFO: Writing fills to: ', outfile)	

    infile = open(infile, 'r')
    outfile = open(outfile, 'w')

    for line in infile:
        line = line.lstrip()

        if line.startswith('net'):
            content = line.split(" ")
            chrom = content[1]

        if line.startswith('fill'):
            content = line.split(" ")
            start = content[1]
            end = str(int(content[1]) + int(content[2]))
            fill_type = content[content.index('type') + 1]

            outfile.write(chrom + "\t" + 
            		start + "\t" + 
            		end + "\t" + 
            		fill_type + "\n")

    infile.close()
    outfile.close()

if __name__ == '__main__':

    # read commandline argumets
    args = commandline()

    checkFillSize(args.input_file, args.output_file)
