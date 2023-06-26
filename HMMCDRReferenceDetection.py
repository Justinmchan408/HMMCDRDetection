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
            prog='HMMCDRReferenceDetection.py',
            description="""Separate candidate CDR containing reads from a bamfile""")
        

        self.parser.add_argument("-m", "--matrix",
                            action = 'store', 
                            nargs='?', 
                            const=True,
                            required=False,
                            default=False,
                            metavar="Runs with Initial CDR Estimate File",
                            help="Runs HMM CDR detection depending if input file CDR Estimate Region file is provided")
        self.parser.add_argument("-a", "--modPosFile",
                            required=True,
                            metavar="bed file containing modified CPG site probabilities",
                            help="The bed file that contains CpG sites and their estimated mod probabilties (File is an output of pb-CpG-Tools/aligned_bam_to_cpg_scores.py)")
        self.parser.add_argument("-b", "--strictBedFile",
                            required=False,
                            metavar="bed file containing estimate CDR Regions",
                            help="bed file containing estimate CDR Regions with chromosome, starting position and ending positon")
        

        self.parser.add_argument("-aa",
                            required=False,
                            type = float,
                            metavar="Probability from next CpG position being in a CDR given current CpG position in CDR",
                            help="Probability from next CpG position being in a CDR given current CpG position in CDR")
        self.parser.add_argument("-ab",
                            required=False,
                            type = float,
                            metavar="Probability from next CpG position not being in a CDR given current CpG position in CDR",
                            help="Probability from next CpG position not being in a CDR given current CpG position in CDR")
        self.parser.add_argument("-bb",
                            required=False,
                            type = float,
                            metavar="Probability from next CpG position not being in a CDR given current CpG position not in CDR",
                            help="Probability from next CpG position not being in a CDR given current CpG position not in CDR")
        self.parser.add_argument("-ba",
                            required=False,
                            type = float,
                            metavar="Probability from next CpG position being in a CDR given current CpG position not in CDR",
                            help="Probability from next CpG position being in a CDR given current CpG position not in CDR")
        
    
        self.parser.add_argument("-ax",
                            required=False,
                            type = float,
                            metavar="Probability of current CpG position is methylated given current CpG position in CDR",
                            help="Probability of current CpG position is methylated given current CpG position in CDR")
        self.parser.add_argument("-ay",
                            required=False,
                            type = float,
                            metavar="Probability of current CpG position is not methylated given current CpG position in CDR",
                            help="Probability of current CpG position is not methylated given current CpG position in CDR")
        self.parser.add_argument("-bx",
                            required=False,
                            type = float,
                            metavar="Probability of current CpG position is methylated given current CpG position not in CDR",
                            help="Probability of current CpG position is methylated given current CpG position not in CDR")
        self.parser.add_argument("-by",
                            required=False,
                            type = float,
                            metavar="Probability of current CpG position is not methylated given current CpG position not in CDR",
                            help="Probability of current CpG position is not methylated given current CpG position not in CDR")
    
        
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class PathEstimate:
    '''
    Creates the initial path estimate from CG positions from reference genome for Viterbi Leanring. Also
    saves other important variable such as CpG positions and their methylation probabilities and 
    current chromosome

    Attributes:
    - cpgSitesHORFile: Subsetted bed file containing CpG sites and probabilties in HOR regions 
    - cdrStrictRegionsFile: Bed file rough estimates of CDR regions

    Parses data in these files to create initial estimates of transition and emission probabilities
    '''
    def __init__ (self, cpgSitesHORFile):
        '''
        Store HOR CpG site probability bed files
        '''
        self.self.cpgSitesHORFile = cpgSitesHORFile
    def getPath (self):
        '''
        Creates the initial path estimate for initial Viterbi Learning

        Addtionally creates global variables such as:
        - self.chrName: Chromosome label/name from bed files
        - self.cpgSitesAndProbs: Dictionary with key as CpG starting postion and value as methylation probability
        '''
        # Creates variables to store CG postions, chromosome names, and estimated path
        self.cpgSitesAndProbs = {}
        self.chrName = ""
        path = ""
        # Opens the CpG HOR methylation probability bed file as input file
        with open(self.cpgSitesHORFile) as inputFile:
            # Loops through every line in the input file
            for line in inputFile:
                # Splits line into columns in a list
                fileColumnSplit = line.strip().split("\t")
                # Saves chromosome name in column 1 of the file
                self.chrName = fileColumnSplit[0]
                # Saves CpG site starting position and methylation probability
                self.cpgSitesAndProbs[int(fileColumnSplit[1])] = float(fileColumnSplit[3])
                # Checks if CpG methylation probability is greater than 50% to store as methylation count
                if float(fileColumnSplit[3]) >= 50:
                    path += "x"
                else:
                    path += "y"
        return path

    def getChrName(self):
        '''
        Returns chromosome name
        '''
        return self.chrName
    
    def getCpGSitesAndProbs(self):
        '''
        Returns dictionary of CpG starting positions as keys and methylation probabilities as values
        '''
        return self.cpgSitesAndProbs



class InitialMatricesEstimate:
    '''
    Creates the initial transition and emission matrices for Viterbi Learning of CDR Regions. Also
    saves other important variable such as CpG positions and their methylation probabilities and 
    current chromosome

    Attributes:
    - cpgSitesHORFile: Subsetted bed file containing CpG sites and probabilties in HOR regions 
    - cdrStrictRegionsFile: Bed file rough estimates of CDR regions

    Parses data in these files to create initial estimates of transition and emission probabilities
    '''
    def __init__ (self, cpgSitesHORFile, cdrStrictRegionsFile):
        '''
        Store CDR region and HOR CpG site probability bed files
        '''
        self.cpgSitesHORFile = cpgSitesHORFile
        self.cdrStrictRegionsFile = cdrStrictRegionsFile

    def getInitialTransitionMatrix (self):
        '''
        Creates the initial transition matrix for CDR Viterbi Learning

        Addtionally creates global variables such as:
        - self.cdrRegions: List of tuples where each tuple represents a CDR region
        - self.chrName: Chromosome label/name from bed files
        - self.cpgSitesAndProbs: Dictionary with key as CpG starting postion and value as methylation probability
        '''

        # Creates variables to store CDR region postions and chromosome names
        self.cdrRegions = []
        self.chrName = ""

        # Opens the Strict CDR File as input file
        with open(self.cdrStrictRegionsFile) as inputFile:
            # Loops through each line in the file
            for line in inputFile:
                # Splits line into columns in a list
                fileColumnSplit = line.strip().split("\t")
                self.chrName = fileColumnSplit[0]
                # Saves CDR Regions from starting and ending positions
                self.cdrRegions.append((int(fileColumnSplit[1]), int(fileColumnSplit[2])))

        # Creates variables to store the initial transition matrix, counts of CpG sites in 
        # CDR/non-CDR positons and CpG sites starting positions along with methylation probabilities
        transitionMatrix = {'AA': 1e-10, 'AB': 1e-10, 'BA': 1e-10, 'BB': 1e-10}
        prevStateCounts = {'A': 2e-10, 'B': 2e-10}
        self.cpgSitesAndProbs = {}

        # Opens the CpG HOR methylation probability bed file as input file
        with open(self.cpgSitesHORFile) as inputFile:
            # Loops through every line in the input file
            for line in inputFile:
                # Splits line into columns in a list
                fileColumnSplit = line.strip().split("\t")
                # Saves CpG site starting position and methylation probability
                self.cpgSitesAndProbs[int(fileColumnSplit[1])] = float(fileColumnSplit[3])

        # Variables to loop through CDR regions based on CpG site position and previous CpG site
        # in a CDR or not in a CDR
        cdrIndex = 0
        prevCPGsiteCDR = None

        # Loops through CpG sites
        for posNum, cpgSitePos in enumerate(self.cpgSitesAndProbs.keys()):
            # Loops until CpG site is in a CDR region or behind it
            while cpgSitePos > self.cdrRegions[cdrIndex][1]:
                cdrIndex += 1
            
            currCDRStartPos = self.cdrRegions[cdrIndex][0]
            currCDREndPos = self.cdrRegions[cdrIndex][1]
            # Checks if CpG site is in a CDR Region
            if cpgSitePos >= currCDRStartPos and cpgSitePos <= currCDREndPos:
                # Checks if this is not the first CpG site in the loop
                if posNum != 0:
                    # Checks if the previous CpG site was in a CDR
                    if prevCPGsiteCDR == "A":
                        transitionMatrix["AA"] += 1
                    # Checks if the previous CpG site was not in a CDR
                    else:
                        transitionMatrix["BA"] += 1
                    prevStateCounts[prevCPGsiteCDR] += 1 # Adds to total count of normalization
                prevCPGsiteCDR = "A"
            # Checks if CpG site is not in a CDR Region
            else:
                # Checks if this is not the first CpG site in the loop
                if posNum != 0:
                    # Checks if the previous CpG site was in a CDR
                    if prevCPGsiteCDR == "A":
                        transitionMatrix["AB"] += 1
                    # Checks if the previous CpG site was not in a CDR
                    else:
                        transitionMatrix["BB"] += 1
                    prevStateCounts[prevCPGsiteCDR] += 1  # Adds to total count of normalization
                prevCPGsiteCDR = "B"

        # Loops through each transition states and counts to normalize
        for transitionStates, transitionStatesCount in transitionMatrix.items():
            transitionMatrix[transitionStates] /= (prevStateCounts.get(transitionStates[0]) + 1e-10)
        
        return transitionMatrix
    
    def getInitialEmissionMatrix (self):
        '''
        Creates the initial emission matrix for CDR Viterbi Learning
        '''

        # Stores emission matrix, # of CpG sites in CDR/not in CDRs, 
        # path/methylation across CpG sites and index of CDR examined for overlap
        emissionMatrix = {'Ax': 1e-10, 'Ay': 1e-10, 'Bx': 1e-10, 'By': 1e-10}
        emissionStateCounts = {'A': 2e-10, 'B': 2e-10}
        path = ""
        cdrIndex = 0

        # Loops through CpG sites
        for cpgSitePos, cpgProb in self.cpgSitesAndProbs.items():
            # Loops until CpG site is in a CDR region or behind it
            while cpgSitePos > self.cdrRegions[cdrIndex][1]:
                cdrIndex += 1
            
            currCDRStartPos = self.cdrRegions[cdrIndex][0]
            currCDREndPos = self.cdrRegions[cdrIndex][1]
            # Checks if CpG site is in a CDR Region
            if cpgSitePos >= currCDRStartPos and cpgSitePos <= currCDREndPos:
                # Checks if CpG methylation probability is greater than 50% to store as methylation count
                if cpgProb >= 50:
                    emissionMatrix["Ax"] += 1
                    path += "x"
                # Checks if CpG methylation probability is less than 50% to store as non-methylation count
                else:
                    emissionMatrix["Ay"] += 1
                    path += "y"
                emissionStateCounts["A"] += 1 # Saves total state count to normalize counts to probabilities
            # Checks if CpG site is not in a CDR Region
            else:
                # Checks if CpG methylation probability is greater than 50% to store as methylation count
                if cpgProb >= 50:
                    emissionMatrix["Bx"] += 1
                    path += "x"
                # Checks if CpG methylation probability is less than 50% to store as non-methylation count
                else:
                    emissionMatrix["By"] += 1
                    path += "y"
                emissionStateCounts["B"] += 1 # Saves total state count to normalize counts to probabilities

        # Loops through emission states and their counts to normalize them to probabilities
        for emission, emissionCount in emissionMatrix.items():
            emissionMatrix[emission] /= emissionStateCounts.get(emission[0])
        
        return path, emissionMatrix
    
    def getChrName(self):
        '''
        Returns chromosome name
        '''
        return self.chrName
    
    def getCpGSitesAndProbs(self):
        '''
        Returns dictionary of CpG starting positions as keys and methylation probabilities as values
        '''
        return self.cpgSitesAndProbs

class ViterbiLearning:
    '''
    Performs Viterbi Learning given an emission path/CpG methylation, transition matrix, and emission
    matrix to return the most probable CDR regions

    Attributes :
    - data: list of various attributes
        - data[0]: string of emission path for CpG sites (X - methylation, Y - not a methylation)
        - data[1]: list of emission states (X - methylation, Y - not a methylation)
        - data[2]: list of transition states (A - CDR Region, B - Not a CDR Region)
        - data[3]: dictionary of transition states and their probabilities
        - data[4]: dictionary of emission states and their probabilities
    '''
    def __init__(self, data):
        '''
        Store path, emission states, transition states, transition matrix and emission matrix as global variables
        '''
        self.path = data[0]
        self.emissionStates = data[1]
        self.transitionStates = data[2]
        self.transitionMatrix = data[3]
        self.emissionMatrix = data[4]

    def parameterEstimation (self, hiddenStates):
        '''
        Performs parameter estimation given an emission path and hidden state path. In other words, recalculates
        transition matrix and emission matrix from the currest emission path and hidden state path estimates

        Attributes:
        - hiddenStates: String of current hidden state path estimate 
        '''

        # Stores new transition and emission matrix along with total counts to get conditional probabilities
        newTransitionMatrix = {transition:1e-10 for transition in self.transitionMatrix.keys()}
        newEmissionMatrix = {emission:1e-10 for emission in self.emissionMatrix.keys()}
        transitionTotalCount = {state:2e-10 for state in self.transitionStates}
        emissionTotalCount = {state:2e-10 for state in self.transitionStates}

        # Loops through emission path
        for currIndex in range(len(self.path)-1):
            # Forms strings that represent current transition state and emission state
            newTransition = hiddenStates[currIndex:currIndex+2]
            newEmission = hiddenStates[currIndex] + self.path[currIndex] 
            transitionTotal = hiddenStates[currIndex]
            emissionTotal = hiddenStates[currIndex]
            
            # Adds counts to transition state and emission states in new matrices
            newTransitionMatrix[newTransition] += 1
            newEmissionMatrix[newEmission] += 1
            transitionTotalCount[transitionTotal] += 1
            emissionTotalCount[emissionTotal] += 1
        
        # Adds count for last emission state in emission matrix 
        newEmission = hiddenStates[-1] + self.path[-1]
        emissionTotal = hiddenStates[-1]
        newEmissionMatrix[newEmission] += 1
        emissionTotalCount[emissionTotal] += 1

        # Loops through transition states and counts to calcuate their probabilities/normalizing
        for transition in self.transitionMatrix.keys():
            newTransitionMatrix[transition] /= transitionTotalCount.get(transition[0])
        
        # Loops through emission states and counts to calcuate their probabilities/normalizing
        for emission in self.emissionMatrix.keys():
            newEmissionMatrix[emission] /= emissionTotalCount.get(emission[0])
        
        return newTransitionMatrix, newEmissionMatrix

    def estimateBestPath (self):
        '''
        Calculates the best/most probable hidden state path from current transition and emission matrices
        using the Viterbi Algorithm. Essentially calculates the most probable path at each CpG site for 
        each ending state (A - CDR, B - Not a CDR)
        '''
        # List of dictionaries that store growing hidden state path at each CpG site and its probability
        viterbiGraph= []
        # Loops through emission path/CpG sites
        for emissionStep in range(len(self.path)):
            # Stores current CpG most probable hidden state path along with their probabilities
            currStepProbs = {}
            # Loops through each hidden state
            for indivState in self.transitionStates:
                # Checks if we are in the first step of the emission path
                if emissionStep == 0:
                    # Saves log values of transition probabilities starting in the beginning
                    currStepProbs[indivState] = np.log(1.0/len(self.transitionStates))
                    emissionStr = indivState + self.path[emissionStep]
                    # Adds log value of emission probability
                    currStepProbs[indivState] += np.log(self.emissionMatrix.get(emissionStr))
                # Checks if we are past the first step of the emission path
                else:
                    bestCurrPath, bestProb = "", float("-inf")
                    # Loops through previous CpG site hidden state paths and their probabilties
                    for prevPath, prevProb in viterbiGraph[emissionStep - 1].items():
                        transmissionStr = prevPath[-1] + indivState
                        emissionStr = indivState + self.path[emissionStep]
                        # Calculates new probabilties of hidden path candidate
                        currTotProb = prevProb + np.log(self.transitionMatrix.get(transmissionStr)) + np.log(self.emissionMatrix.get(emissionStr))
                        # Checks if hidden path candidate probability is the highest
                        if currTotProb > bestProb:
                            bestProb = currTotProb
                            bestCurrPath = prevPath + indivState
                    # Saves most probable hidden state path and its probability of current CpG
                    currStepProbs[bestCurrPath] = bestProb
            viterbiGraph.append(currStepProbs)
        # Returns most hidden path once the last two candidates are generated
        return max(viterbiGraph[len(viterbiGraph)-1], key=viterbiGraph[len(viterbiGraph)-1].get)



    def performViterbiLearning(self):
        '''
        Perform Viterbi Learning until convergence of transition and emission matrices. Involves a loop
        to use emission path and transition & emission matrices to estimated best hidden state path and
        then emission path and the estimated hidden path to calculate new transition & emission matrices
        '''

        # Creates a set of transition and emission probabilties to compare if they are changing
        prevTransitionSum = set(self.transitionMatrix.values())
        prevEmissionSum = set(self.emissionMatrix.values())
        currTransitionSum = set()
        currEmissionSum = set()
        
        # Loops until the transition and emission matrices stay the same/converge
        while prevTransitionSum != currTransitionSum or prevEmissionSum != currEmissionSum:
            prevTransitionSum = set(self.transitionMatrix.values())
            prevEmissionSum = set(self.emissionMatrix.values())

            hiddenStates = self.estimateBestPath() # Viterbi Algorithm
            self.transitionMatrix, self.emissionMatrix = self.parameterEstimation(hiddenStates) # Estimate HMM Parameters

            currTransitionSum = set(self.transitionMatrix.values())
            currEmissionSum = set(self.emissionMatrix.values())

        return self.estimateBestPath()

    def generateBedFile(self, chr, estimatedStates, cpgSitesAndProbs):
        '''
        Creates bed file that contains CDR regions found through Viterbi Learning

        Attributes:
        - chr: string of chromosome name/label
        - estimatedStates: string offinal hidden state path found in Viterbi Learning
        - cpgSitesAndProbs: dictionary of CpG starting positions as keys and their methylation probabilties as values
        '''

        newCDRRegions = [] # List to store starting and stopping positions of Viterbi Learning CDRs

        # Stores if previous and current CpG sites were in a CDR or not in a CDR
        prevState = None
        prevPos = None
        hitStart = False

        # Loops through hidden state path and starting positions of CpG sites
        for stateEst, cpgPos in zip(estimatedStates, cpgSitesAndProbs.keys()):
            currState = stateEst
            # Checks for starting position for a CDR
            if prevState == "B" and currState == "A":
                newCDRRegions.append([chr, cpgPos, -1])
                hitStart = True
            # Checks for and ending position for a CDR
            elif prevState == "A" and currState == "B":
                # Checks if a start to a CDR has not been encountered
                if hitStart == False:
                    newCDRRegions.append([chr, list(cpgSitesAndProbs.keys())[0], prevPos])
                # Checks if a start to a CDR has been encountered
                else:
                    newCDRRegions[len(newCDRRegions)-1][2] = prevPos
            prevState = stateEst
            prevPos = cpgPos
        
        # Checks if there was no proper ending to a CDR
        if hitStart and newCDRRegions[len(newCDRRegions)-1][2] == -1:
            newCDRRegions[len(newCDRRegions)-1][2] = list(cpgSitesAndProbs.keys())[len(cpgSitesAndProbs)-1]
        
        # Loops through CDR starting and stopping positions to format an output bed file of these CDR regiosn
        for cdr in newCDRRegions:
            print(cdr[0] + "\t" + str(cdr[1]) + "\t" + str(cdr[2]))






def main(options=None):
    '''
    Initializes the first transition and emission matrices and then runs Viterbi Learning to estimate CDR Regions
    '''
    # Creates command line obejct
    thisCommandLine = None
    # Checks if options/parameters were entered
    if options is None:
        thisCommandLine = CommandLine()
    else:
        thisCommandLine = CommandLine(options)
    
    # Runs if initial matrix file is provided
    if thisCommandLine.args.matrix:
        # Creates InitialMatricesEstimate object to get initial transition and emission matrix
        matrixEstimator = InitialMatricesEstimate(thisCommandLine.args.modPosFile, thisCommandLine.args.strictBedFile)
        initialTransitionMatrix = matrixEstimator.getInitialTransitionMatrix()
        path, initialEmissionMatrix = matrixEstimator.getInitialEmissionMatrix()
        chrName = matrixEstimator.getChrName()
        cpgSitesAndProbs = matrixEstimator.getCpGSitesAndProbs()
    else:
        initialTransitionMatrix = {'AA': thisCommandLine.args.aa, 'AB': thisCommandLine.args.ab, 'BA': thisCommandLine.args.ba, 'BB': thisCommandLine.args.bb}
        initialEmissionMatrix = {'Ax': thisCommandLine.args.ax, 'Ay': thisCommandLine.args.ay, 'Bx': thisCommandLine.args.bx, 'By': thisCommandLine.args.by}
        # Creates pathEstimator object to get initial transition and emission matrix
        pathEstimator = PathEstimate(thisCommandLine.args.modPosFile)
        path = pathEstimator.getPath()
        chrName = pathEstimator.getChrName()
        cpgSitesAndProbs = pathEstimator.getCpGSitesAndProbs()

    # States and emission of HMM Model
    states = ["A", "B"] # A - CDR, B - Not a CDR
    emissions = ["x", "y"] # X - Methylation, Y - No Methylation

    # Creates a list for input for the ViterbiLearning object to output optimal CDR regions in a bed file
    viterbiData = [path, emissions, states, initialTransitionMatrix, initialEmissionMatrix]
    vitLearn = ViterbiLearning(viterbiData)
    estimatedStates = vitLearn.performViterbiLearning()
    vitLearn.generateBedFile(chrName, estimatedStates, cpgSitesAndProbs)


if __name__ == '__main__':
    main()