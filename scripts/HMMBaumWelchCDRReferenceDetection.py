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
            prog='HMMBaumWelchCDRReferenceDetection.py',
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
        self.parser.add_argument("-b", "--initialCDRBedFile",
                            required=False,
                            metavar="bed file containing estimate CDR Regions",
                            help="bed file containing estimate CDR Regions with chromosome, starting position and ending positon")
        self.parser.add_argument("-c", "--modCutoff",
                            required=False,
                            type = float,
                            default=50.0,
                            metavar="Methylation percentage cutoff for modification",
                            help="A float that represents a methylation percentage cutoff. Mod site metyhlation prob >= modCutoff(methylated), Mod site metyhlation prob < modCutoff (not methylated) ")
        self.parser.add_argument("-o", "--matrixOutputFile",
                            required=True,
                            metavar="txt file containing initial matrices and final matrices",
                            help="A file that is written initial and final transition and emission matrices of the Baum Welch algorithm")
        

        self.parser.add_argument("-aa",
                            required=False,
                            type = float,
                            default=98.6981,
                            metavar="Probability from next CpG position being in a CDR given current CpG position in CDR",
                            help="Probability from next CpG position being in a CDR given current CpG position in CDR")
        self.parser.add_argument("-ab",
                            required=False,
                            type = float,
                            default=1.3019,
                            metavar="Probability from next CpG position not being in a CDR given current CpG position in CDR",
                            help="Probability from next CpG position not being in a CDR given current CpG position in CDR")
        self.parser.add_argument("-bb",
                            required=False,
                            type = float,
                            default=99.9368,
                            metavar="Probability from next CpG position not being in a CDR given current CpG position not in CDR",
                            help="Probability from next CpG position not being in a CDR given current CpG position not in CDR")
        self.parser.add_argument("-ba",
                            required=False,
                            type = float,
                            default=0.0632,
                            metavar="Probability from next CpG position being in a CDR given current CpG position not in CDR",
                            help="Probability from next CpG position being in a CDR given current CpG position not in CDR")
        
    
        self.parser.add_argument("-ax",
                            required=False,
                            type = float,
                            default=11.1472,
                            metavar="Probability of current CpG position is methylated given current CpG position in CDR",
                            help="Probability of current CpG position is methylated given current CpG position in CDR")
        self.parser.add_argument("-ay",
                            required=False,
                            type = float,
                            default=88.8528,
                            metavar="Probability of current CpG position is not methylated given current CpG position in CDR",
                            help="Probability of current CpG position is not methylated given current CpG position in CDR")
        self.parser.add_argument("-bx",
                            required=False,
                            type = float,
                            default=81.0364,
                            metavar="Probability of current CpG position is methylated given current CpG position not in CDR",
                            help="Probability of current CpG position is methylated given current CpG position not in CDR")
        self.parser.add_argument("-by",
                            required=False,
                            type = float,
                            default=18.9636,
                            metavar="Probability of current CpG position is not methylated given current CpG position not in CDR",
                            help="Probability of current CpG position is not methylated given current CpG position not in CDR")
    
        
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class HMMCDRInitializer:
    '''
    Creates initial emission path estimate, initial transition and emission matrix to begin Baum-Welch learning for CDRs.
    The class is made to set up for processing the possible modification sites for analysis for the Baum-Welch Algorithm and
    return those variables.

    Attributes:
    - cpgSitesHORFile: Subsetted bed file containing modification sites and probabilties in HOR regions 
    - (optional) cdrInitialRegionsFile: Bed file with initial estimates of CDR regions

    Parses data in these files to create initial estimates of transition and emission probabilities
    '''
    def __init__ (self, modSiteFile):
        '''
        Saves modification methylation file
        '''
        self.modSiteFile = modSiteFile

    def getPath(self, modCutoff = 50):
        '''
        Generates an emissions path from methylation file

        Loops through each possible modification site and saves each site as methylated (0) or not (1) based
        on methylation file

        Assumptions:
        - Methylation file has at least four columns
            - (1) Chromosome Name
            - (2) Left Position of Site
            - (3) Right Position of Site
            - (4) Methylation percentage (Range: 0-100)
        
        Output:
        - Emissions path
        - A list is lists where each sublist contains chromosome name,
          left and right position, and methylation of each site
        '''
        # Stores emissions path and modification site data
        self.path = []
        self.modPositions = []

        # Open the mod input file
        with open(self.modSiteFile) as modInputFile:
            # Loops through each line in input file
            for line in modInputFile:
                # Attempts to run but will skip over sites with n/a, nan in methylation percentage
                try:
                    # Saves chromosome name, left and right position, and methylation percentage
                    splitLine = line.split("\t")
                    chrName = splitLine[0]
                    firstCoord = int(splitLine[1])
                    lastCoord = int(splitLine[2])
                    modProb = float(splitLine[3])
                    # Checks if mod site has methylation higher than threshold
                    if modProb >= modCutoff:
                        self.path.append(0)
                        self.modPositions.append([chrName, firstCoord, lastCoord, modProb, 0])
                    # Methylation is lower than threshold
                    else:
                        self.path.append(1)
                        self.modPositions.append([chrName, firstCoord, lastCoord, modProb, 1])
                except:
                    continue

        self.path = np.array(self.path) # Saves emissions path as array
        return self.path, self.modPositions

    def getStatesAndEmissions(self):
        '''
        Generates arrays of possible transitions and emissions
        
        Output:
        - Possible states array (0 - CDR, 1 - Not a CDR)
        - Possible emissions array (0 - Methylated, 1 - Not methylated)
        '''
        self.states = np.array([0, 1])
        self.emissions = np.array([0, 1])
        return self.states, self.emissions
    
    def getMatricesWithCDRFile(self, cdrEstimateFile):
        '''
        Generates an initial transition and emission matrices as arrays from bed file

        Loops through initial CDR estimate bed file to find initial guesses
        for transition and emission probabilities

        Assumptions:
        - Initial CDR Estimate file has at least three columns
            - (1) Chromosome Name
            - (2) Left Position of CDR Region
            - (3) Right Position of CDR Region
        
        Output:
        - Transition Matrix Array
        - Emissions Matrix Array
        '''
        cdrEstimateRegions = [] # Save CDR Initial Estimates
        # Open the mod input file
        with open(cdrEstimateFile) as cdrFile:
            # Loops through each line in the file
            for line in cdrFile:
                # Create variables for each columns in the file
                splitLine = line.split("\t")
                chrName = splitLine[0]
                firstCdrCoord = int(splitLine[1])
                lastCdrCoord = int(splitLine[2])
                # Saves CDR Regions
                cdrEstimateRegions.append([chrName, firstCdrCoord, lastCdrCoord])
        
        # Creates initial transition and emission matrices as arrays
        self.transitionMatrix = np.zeros((len(self.states), len(self.states)))
        self.emissionMatrix = np.zeros((len(self.states), len(self.emissions)))
        
        cdrRegionIndex = 0
        prevModState = None

        # Loops through possible mod sites
        for chrName, firstModCoord, secondModCoord, modProb, emission in self.modPositions:
            cdrLastPos = cdrEstimateRegions[cdrRegionIndex][2]
            # Loops through CDR regions until last position of CDR region comes after mod site position
            while cdrLastPos <= firstModCoord:
                cdrRegionIndex += 1
                cdrLastPos = cdrEstimateRegions[cdrRegionIndex][2]
            
            cdrBegPos = cdrEstimateRegions[cdrRegionIndex][1]
            withinCDR = 1

            # Checks if possible mod position is in CDR
            if cdrBegPos <= firstModCoord:
                withinCDR = 0
            
            # Saves counts for emissions matrix
            self.emissionMatrix[withinCDR, emission] += 1

            # Saves counts for transition matrix
            if prevModState is not None:
                self.transitionMatrix[prevModState, withinCDR] += 1
            prevModState = withinCDR # Updates previous mod states (0 - CDR, 1 - Not CDR)

        # Normalizes transition and emissions matrices into probabilities
        transitionStateSum = np.sum(self.transitionMatrix, axis=1)
        self.transitionMatrix = np.transpose(np.divide(np.transpose(self.transitionMatrix), transitionStateSum))
        emissionSum = np.sum(self.emissionMatrix, axis=1)
        self.emissionMatrix = np.transpose(np.divide(np.transpose(self.emissionMatrix), emissionSum))

        return self.transitionMatrix, self.emissionMatrix


    def getMatricesFromCommandLine(self, aa, ab, ba, bb, ax, ay, bx, by):
        '''
        Generates an initial transition and emission matrices as arrays

        Assumptions:
        - All probabilties should range from 0-100
            - aa + ab = 100
            - ba + bb = 100
            - ax + ay = 100
            - bx + by = 100
        
        Output:
        - Transition Matrix Array
        - Emissions Matrix Array
        '''
        self.transitionMatrix = np.array([[aa/100, ab/100],[ba/100, bb/100]])
        self.emissionMatrix = np.array([[ax/100, ay/100],[bx/100, by/100]])
        return self.transitionMatrix, self.emissionMatrix

class HMMBaumWelchCDRDetection:
    '''
    Performs Baum Welch Learning to give probabilities if possible mod site in a CDR

    Attributes:
    - path: array emission path of possible mod sites (0 - Methylated, 1 - Not methylated)
    - transition matrix: array of transition probabilities between CDR/Non CDR regions
    - emission matrix: array of emission probabilities of methylation/no methylation given in a particular region
    - transitions: transition states (0 - CDR, 1 - CDR)
    - emissions: emission states (0 - Methylated, 1 - Not methylated)

    Alternates between M-step and E-step of algorithm to determine probability if site in CDR
    '''
    def updateModCDRProbabilties(self, path, emissions, states, realTrans, realEm):
        '''
        Calculates pi star and pi star star matrix from Baum Welch Learning

        Performs Forward and Backward Algorithm and uses matrices storing values
        to be able to calculate forward and backwards algorithm.
        
        Output:
        - Pi Star Matrix: Calculates probability of reference mod sites in CDR or not
        - Pi Star Star (Double) Matrix: Calculate probability of each 
          transition from one site to the next
        '''
        # Creates arrays for forward and backward algorithm sums
        forwardAlgo = np.zeros((len(states), len(path)))
        backwardsAlgo = forwardAlgo.copy()
        
        # Initializes initial array for forward and backward algorithm
        forwardAlgo[:,0] = np.log(1/len(states)) + np.log(realEm[:,path[0]])

        # Loops through each site to perform Forward and Backwards Algorithm
        for i in range(1, len(path)):
            # Forward Algorithm
            a = np.transpose(forwardAlgo[:,i-1] + np.transpose(np.log(realTrans) + np.log(realEm[:,path[i]])))
            b = np.max(a, axis=0) # Log-sum-exp trick for underflow
            forwardAlgo[:,i] = b + np.log(np.sum(np.exp(a-b), axis=0))

            # Backwards Algorithm
            a = backwardsAlgo[:,-i] + np.log(realTrans) + np.log(realEm[:,path[-i]])
            b = np.max(a, axis=1) # Log-sum-exp trick for underflow
            backwardsAlgo[:,-i-1] = b + np.log(np.sum(np.exp(np.transpose(np.transpose(a)-b)), axis=1))


        # Perform calculations for pi star and pi double star matrices
        forwardBackward = forwardAlgo + backwardsAlgo
        sinkLogSumMax = np.max(forwardAlgo[:,-1])
        sinkLogSum = sinkLogSumMax + np.log(np.sum(np.exp(forwardAlgo[:,-1] - sinkLogSumMax)))

        # Pi Star
        piStar = np.exp(forwardBackward-sinkLogSum)

        piDoubleStar = np.zeros((len(states) * len(states), len(path)-1))
        # Loops through each transition between sites to calculate Pi Double Star
        for j in range(0, len(path)-1):
            piDoubleStar[:,j] = (np.add.outer(forwardAlgo[:,j], backwardsAlgo[:,j+1]) + np.log(realTrans) + np.log(realEm[:,path[j+1]])).flatten() - sinkLogSum 
        piDoubleStar = np.exp(piDoubleStar)

        return piStar, piDoubleStar

    def updateMatrices(self, path, states, emissions, piStar, piDoubleStar):
        '''
        Calculates the updated transition and emission matrices for Baum Welch Algorithm

        Nomalizes the Pi Star and Pi Double Star matrices to generate updated matrices
        
        Output:
        - Transition Matrix: Calculates updated transition matrix from Pi Double Star Matrix
        - Emission Matrix: Calculates updated emission matrix from Pi Star Matrix
        '''
        # Generate updated Transition Matrix 
        transitionStateSum = np.sum(piDoubleStar, axis=1).reshape(2,2)
        transitionInitSum = np.sum(transitionStateSum, axis=1)
        transitionMatrix = np.transpose( np.transpose(transitionStateSum) / transitionInitSum)

        # Generate updated Emission Matrix
        emissionStateSum = np.sum(piStar, axis=1)
        emissionMatrix = np.zeros((len(states), len(emissions)))
        for k in range(0, len(path)):
            emissionMatrix[:,path[k]] += piStar[:,k]
        emissionMatrix = np.transpose(np.transpose(emissionMatrix) / emissionStateSum)

        return transitionMatrix, emissionMatrix

    def outputBedFile(self, piStar, modPositions):
        '''
        Generates an output bed file to visualize regions in IGV
        
        Output:
        - Stdout: Prints CDR Predictions from Baum Welch Algorithm to bed file
        '''

        # Begins parsing and storing first mod in condensed bed file
        chrName, firstModCoord, secondModCoord, modProb, emission = modPositions[0]
        cdrProb = round(piStar[0][0] * 100, 2)
        # List containing CDR Regions
        condensedModPos = [[chrName, firstModCoord, secondModCoord, cdrProb]]

        # Loops through every mod
        for modIndex in range(1, piStar.shape[1]):
            chrName, firstModCoord, secondModCoord, modProb, emission = modPositions[modIndex]
            # Check if the methylation percentage is the same to the second decimal place
            if round(piStar[0][modIndex] * 100, 2) == condensedModPos[-1][3]:
                condensedModPos[-1][2] = secondModCoord # Builds on from past mod if methylation is the similar
            else:
                # Creates new regions if methylation differs from previous mod
                condensedModPos.append([chrName, firstModCoord, secondModCoord, round(piStar[0][modIndex] * 100, 2)])
        
        print("track name=\"HMM Baum Welch CDR\" description=\"HMM Baum Welch CDR Regions\" itemRgb=\"On\"")
        
        # Loops through the condensed mod regions
        for modIndex in range(len(condensedModPos)):
            chrName, firstModCoord, secondModCoord, cdrProb = condensedModPos[modIndex]
            # Increases bounds of regions past mod sites
            
            # Will expand first/beginning coordinate of CDR if not first condensed mod region
            if modIndex > 0:
                firstModCoord = round((firstModCoord + condensedModPos[modIndex-1][2])/2)
            # Will expand last/end coordinate of CDR if not last condensed mod region
            if modIndex < len(condensedModPos) - 1:
                secondModCoord = round((secondModCoord + condensedModPos[modIndex+1][1])/2)
            
            # Generates RGB value based on CDR prediction
            red, green, blue = round(255 - 255 * (cdrProb/100)), round(255 - 255 * (cdrProb/100)), 255
            color = str(red) + "," + str(green) + "," + str(blue)
            print(*[chrName, firstModCoord, secondModCoord, "CDR Region", 0, ".", firstModCoord, secondModCoord, color,  cdrProb], sep='\t')
        

def main(options=None):
    '''
    Initializes the first transition and emission matrices and then runs Baum Welch Algorithm to estimate CDR Regions
    '''
    # Creates command line obejct
    thisCommandLine = None
    # Checks if options/parameters were entered
    if options is None:
        thisCommandLine = CommandLine()
    else:
        thisCommandLine = CommandLine(options)

    # Creates HMMCDRInitializer object to create initial variables needed to begin Baum-Welch Algorithm
    hmmInit = HMMCDRInitializer(thisCommandLine.args.modPosFile)
    # Returns emission path and sites of possible methylation
    path, modPositions = hmmInit.getPath(thisCommandLine.args.modCutoff)
    # Returns numpy array of states: [0, 1](CDR, not CDR) and [0, 1](methylated, not methylated)
    states, emissions = hmmInit.getStatesAndEmissions()

    # Runs if initial matrix file is provided
    # Creates transition and emission matrices
    if thisCommandLine.args.matrix:
        transitionMatrix, emissionMatrix = hmmInit.getMatricesWithCDRFile(thisCommandLine.args.initialCDRBedFile)
    else:
        transitionMatrix, emissionMatrix = hmmInit.getMatricesFromCommandLine(thisCommandLine.args.aa, thisCommandLine.args.ab, thisCommandLine.args.ba, thisCommandLine.args.bb, thisCommandLine.args.ax, thisCommandLine.args.ay, thisCommandLine.args.bx, thisCommandLine.args.by)

    # Stores previous transition and emission matrices from previous iteration
    prevTransitionMatrix, prevEmissionMatrix = np.zeros(transitionMatrix.shape), np.zeros(emissionMatrix.shape)
    
     # Creates HMMBaumWelchCDRDetection object to perform the Baum-Welch Algorithm to detect CDRs
    hmmCDRDetector = HMMBaumWelchCDRDetection()

    # Write initial matrices lines
    matricesWrite = ["Cutoff: " + str(thisCommandLine.args.modCutoff) + "\n"]
    matricesWrite = matricesWrite + ["Initial Transition Matrix: \n", str(transitionMatrix[0,0]) + "\t", str(transitionMatrix[0,1]) + "\n", str(transitionMatrix[1,0]) + "\t", str(transitionMatrix[1,1]) + "\n"]
    matricesWrite = matricesWrite + ["Initial Emission Matrix: \n", str(emissionMatrix[0,0]) + "\t", str(emissionMatrix[0,1]) + "\n", str(emissionMatrix[1,0]) + "\t", str(emissionMatrix[1,1]) + "\n"]

    # Loops until the difference between all transitions and emissions is less than 0.0001
    while np.max(transitionMatrix-prevTransitionMatrix) > 0.0001 or np.max(emissionMatrix-prevEmissionMatrix) > 0.0001:
        # Sets previous transition and emission matrices
        prevTransitionMatrix = transitionMatrix
        prevEmissionMatrix = emissionMatrix

        # Runs E step: Generates pi star and pi star star matrices using forward and backward algorithms
        piStar, piDoubleStar = hmmCDRDetector.updateModCDRProbabilties(path, emissions, states, transitionMatrix, emissionMatrix)
        # Runs M step: Generates updated transition and emission matrices
        transitionMatrix, emissionMatrix = hmmCDRDetector.updateMatrices(path, states, emissions, piStar, piDoubleStar)

    # Write final matrices lines
    matricesWrite = matricesWrite + ["Final Transition Matrix: \n", str(transitionMatrix[0,0]) + "\t", str(transitionMatrix[0,1]) + "\n", str(transitionMatrix[1,0]) + "\t", str(transitionMatrix[1,1]) + "\n"]
    matricesWrite = matricesWrite + ["Final Emission Matrix: \n", str(emissionMatrix[0,0]) + "\t", str(emissionMatrix[0,1]) + "\n", str(emissionMatrix[1,0]) + "\t", str(emissionMatrix[1,1]) + "\n\n"]

    # Write to matrix file
    matrixFile = open(thisCommandLine.args.matrixOutputFile, "a")
    matrixFile.writelines(matricesWrite)
    matrixFile.close()
    
    # Runs E step one more time to use pi star to determine probability of site to be in CDR or not
    piStar, piDoubleStar = hmmCDRDetector.updateModCDRProbabilties(path, emissions, states, transitionMatrix, emissionMatrix)

    # Creates bed file to visualize in IGV
    hmmCDRDetector.outputBedFile(piStar, modPositions)
    

if __name__ == '__main__':
    main()
