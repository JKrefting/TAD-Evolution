''' 
This script parses a chain alignment file 
(https://genome.ucsc.edu/goldenpath/help/chain.html) and writes only the whole 
chains of the target (reference) species to output file in BED format.
'''

import argparse

def chain2bed(infile, outfile):

    print('INFO: Chain file : ', infile)
    print('INFO: Writing chains to: ', outfile)

    with open(outfile, 'w') as outHandle:
      with open(infile) as inHandle:
        for line in inHandle:
          
          if line.startswith("chain"):
            
            fields = line.strip().split()
            chrom = fields[2]
            strand = fields[4]
            start = fields[5]
            end = fields[6]
            score = fields[1]
            name = fields[12]
            
            outLine = "\t".join([chrom, start, end, name, score, strand])
            outHandle.write(outLine + "\n")

def commandline():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Path to input file in  chain format. See https://genome.ucsc.edu/goldenpath/help/chain.html")
    parser.add_argument("-o", "--output_file", type=str, required=False, default="outfile.chain.bed", help="Path to output file.")
    return parser.parse_args()


if __name__ == '__main__':

    # read commandline argumets
    args = commandline()

    chain2bed(args.input_file, args.output_file)

  
