import numpy as np
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
            prog='ONTSimplifier.py',
            description="""Condenses modbam2bed output by averaging CpG sites on different strands""")
        
        self.parser.add_argument("-a", "--modPosFile",
                            required=True,
                            metavar="bed file containing modified CPG site probabilities",
                            help="The bed file that contains CpG sites and their estimated mod probabilties (File is an output of epi2me-labs/modbam2bed)")
    
        
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class OntModSimplier:
    def __init__(self, cpgSitesHORFile):
        '''
        Store CDR region and HOR CpG site probability bed files
        '''
        self.cpgSitesHORFile = cpgSitesHORFile

    def simplifyFile(self):
        '''
        Concats mods on different strands for ONT
        '''
        with open(self.cpgSitesHORFile) as inputFile: # Open the input file
            lineC = 0 # stores line count
            placeHolder = [] # holds data to merge mod data
            # Loops through each line in the code
            for line in inputFile:
                splitLine = line.split("\t")
                # Checks if currently on even line in the text
                if lineC % 2 == 0:
                    placeHolder = [splitLine[0], splitLine[1], splitLine[2], splitLine[10]]
                else:
                    # Filters out nan as much as possible by taking the average of the mods
                    if placeHolder[3] == "nan" and splitLine[10] == "nan":
                        avgMod = "nan"
                    elif placeHolder[3] == "nan":
                        avgMod = float(splitLine[10])
                    elif splitLine[10] == "nan":
                        avgMod = float(placeHolder[3])
                    else:
                        avgMod = (float(placeHolder[3]) + float(splitLine[10]))/2
                    placeHolder[3] = avgMod
                    print(*placeHolder, sep='\t')

                lineC += 1

def main(options=None):
    '''
    Cleans up ONT mod file for CpG sites on different/opposite strands
    '''

    # Creates command line obejct
    thisCommandLine = None
    # Checks if options/parameters were entered
    if options is None:
        thisCommandLine = CommandLine()
    else:
        thisCommandLine = CommandLine(options)
    
    # Creates OntModSimplier which takes average methylation score on different strands
    modSimple = OntModSimplier(thisCommandLine.args.modPosFile)
    modSimple.simplifyFile()

if __name__ == '__main__':
    main()