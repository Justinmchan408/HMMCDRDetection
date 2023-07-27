# HMM_CDR_detect
Detecting CDRs in HiFi and ONT long reads with methylation tags using Baum Welch Algorithm

```
|__ scripts
    |__ HMMCDRReferenceDetection.py returns CDR regions in HOR centromere regions using mod methylation probabilities with Viterbi Learning
    |__ HMMCDRBaumWelchReferenceDetection.py returns CDR regions in HOR centromere regions using methylation probabilities with the Baum Welch Algorithm
```

## How to run HMMBaumWelchCDRReferenceDetection.py

Python libraries needed:
```
- argparse
- numpy
```

Example Command w/ Initial CDR Region Estimates:

```
python3 HMMBaumWelchCDRReferenceDetection.py -m \
-a mod_pos_file.bed -b inital_cdr_estimate_file.bed \
> HMMCDRBaumWelchOutput.bed
```

Example Command w/ Initial Transition and Emission Probabilities:
```
python3 HMMBaumWelchCDRReferenceDetection.py -a mod_pos_file.bed \
-aa 99.00 -ab 1.00 -ba 1.00 -bb 99.00 \
-ax 20.00 -ay 80.00 -bx 80.00 -by 20.00 \
> HMMCDRBaumWelchOutput.bed
```

Inputs:
- `-a` (required) bed file containing methylation probabilties. A tab delimited file (.tsv) with chromosome name as first column, start position of each site as second column, and methylation probabilty as a percent in the fourth column
- `-b` (optional) bed file containing inital CDR region estimates in centromere. A tab delimited column (.tsv) withg chromosome as first column, start position of CDR region in second column and end position of CDR region of last column
- `-m` (default = False, optional) Determines method to use inout file for CDR regions initial transition and emission matrices or entered values. Default is using entered values for the matrices
- `-aa` (optional, float) Percentage from next CpG position being in a CDR given current CpG position in CDR
- `-ab` (optional, float) Percentage from next CpG position not being in a CDR given current CpG position in CDR
    - -aa and -ab should add up to 100.00
- `-ba` (optional, float) Percentage from next CpG position not being in a CDR given current CpG position not in CDR
- `-bb` (optional, float) Percentage from next CpG position being in a CDR given current CpG position not in CDR
    - -ba and -bb should add up to 100.00
- `-ax` (optional, float) Percentage of current CpG position is methylated given current CpG position in CDR
- `-ay` (optional, float) Percentage of current CpG position is not methylated given current CpG position in CDR
    - -ax and -ay should add up to 100.00
- `-bx` (optional, float) Percentage of current CpG position is methylated given current CpG position not in CDR
- `-by` (optional, float) Percentage of current CpG position is not methylated given current CpG position not in CDR
    - -bx and -by should add up to 100.00


Description of method:

- Begins with an estimated path of methylation of CG sites (or any other methylated sites) within the centromere regions
    - The algorithm also begins with initial transition matrices and emission matrices calculated from a file or inputted      on the command line
- Uses an Baum-Welch EM algorithm using a first order HMM to determine CDRs
    -   E step: Estimate pi matrix (Probabities of each site in a CDR/Not in CDR) and double pi matrix (Probabities of each edge when moving from previous mod site to next mod site)
    -   M step: Resestimate transition and emission matrices by normalizing over pi and double pi matrix
- Continues running until algorithm convergences (When matrices probabilties are changining by less than 0.0001
- Return bed file containing chromosome name, start position and end position of CDR estimates, their probabilities of being a CDR 
