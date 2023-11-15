import argparse

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        self.parser = argparse.ArgumentParser(
            prog='HMMCDRReferenceDetection.py',
            description="""Separate candidate CDR containing reads from a bamfile""")
        

        self.parser.add_argument("-t", "--threshold",
                            required=True,
                            metavar="Strict CDR likelihood percentage threshold to be considered a valid CDR",
                            default=50,
                            help="The percentage threshold to strictly define CDR regions from the Baum Welch Algorithm")
        self.parser.add_argument("-b", "--baumWelchCDRFile",
                            required=True,
                            metavar="Bed file containing estimated likelihood CDR Regions from HMMBaumWelchCDRReferenceDetection.py",
                            help="Bed file containing estimate CDR Regions with chromosome, starting position, ending positon, and CDR likelihood percentages")
     
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)
    
class BaumWelchStrictCDRRegions:
    '''
    Creates strict CDR estimates based on Baum Welch percentage estimates of CDR methylation data.
    Prints starts and ends of estimated CDR Regions given the likelihood of methylation data.

    Attributes :
    - baumWelchFile: bed file containing genomic regions their CDR percentage estimates
    - threshold: percentage threshold for a region to be considered a CDR
    '''
    def __init__(self, baumWelchFile, threshold = 90):
        '''
        Stores Baum Welch CDR Estimate File and strict CDR threshold as global variables
        '''
        self.baumWelchFile = baumWelchFile
        self.threshold = threshold
    
    def createStrictCDRRegions(self):
        '''
        Provides a bed file containing strict CDR estimates from Baum Welch estimates given a threshold

        Loops through each line finding when a CDR region is above and below the threshold to determine
        the start and end of CDR.
        '''
        cdrRegion = False
        leftCDRcoor, maxPerc = 0, 0
        # Opens bed file containing Baum Welch CDR Estimates
        with open(self.baumWelchFile) as inputFile:
            # Loops through each line in the file
            for line in inputFile:
                fileColumnSplit = line.strip().split("\t")
                # Check if its a valid CDR line
                if len(fileColumnSplit) == 10:
                    chr, leftCoor, rightCoor, cdrPerc = fileColumnSplit[0], fileColumnSplit[1], fileColumnSplit[2], float(fileColumnSplit[9])

                    # Checks if currently in a CDR region
                    if cdrRegion:
                        if cdrPerc > maxPerc:
                            maxPerc = cdrPerc
                        # Prints once it hits a CDR region that is smaller than threshold
                        if cdrPerc < self.threshold:
                            print(chr + "\t" + leftCDRcoor + "\t" + leftCoor + '\t' + str(maxPerc) + "\n")
                            cdrRegion = False
                    
                    # Handles not in CDRs
                    else:
                        # Looks for start of CDR from Baum Welch Estimates above threshold
                        if cdrPerc >= self.threshold:
                            cdrRegion = True
                            leftCDRcoor = leftCoor
                            maxPerc = cdrPerc


def main(options=None):
    '''
    Creates BaumWelchStrictCDRRegions object to provide hard estimates for CDRs given a threshold
    '''
    # Creates command line obejct
    thisCommandLine = None
    # Checks if options/parameters were entered
    if options is None:
        thisCommandLine = CommandLine()
    else:
        thisCommandLine = CommandLine(options)

    # Creates a BaumWelchStrictCDRRegions object to provide hard CDR estimates
    aggregator = BaumWelchStrictCDRRegions(thisCommandLine.args.baumWelchCDRFile, thisCommandLine.args.threshold)
    aggregator.createStrictCDRRegions()


if __name__ == '__main__':
    main()

    

                        