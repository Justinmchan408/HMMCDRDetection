# HMM_CDR_detect
Detecting CDRs in HiFi and ONT long reads with methylation tags using Viterbi Learning

```
|__ scripts
    |__ HMMCDRReferenceDetection.py returns CDR regions in HOR centromere regions using mod methylation probabilities with Viterbi Learning
    |__ HMMCDRBaumWelchReferenceDetection.py returns CDR regions in HOR centromere regions using methylation probabilities with Baum Welch Learning
```

## How to run HMMCDRReferenceDetection.py

Python libraries needed:
```
- argparse
- numpy
```

Command w/ Initial CDR Region Estimates:
```
python3 HMMCDRReferenceDetection.py -a mod_pos_bed_file -b inital_cdr_estimate_bed_file > outfile.bed
```

Command w/ Initial Transition and Emission Probabilities:
```
python3 HMMCDRReferenceDetection.py -a mod_pos_bed_file -m -aa AA_prob -ab AB_prob -ba BA_prob -bb BB_prob -ax AX_emission -ay AY_emission -bx BX_emission -by BY_emission > outfile.bed
```

Inputs:
- `-a` (required) bed file containing methylation probabilties. A tab delimited file (.tsv) with chromosome name as first column, start position of each site as second column, and methylation probabilty as a percent in the fourth column
- `-b` (optional) bed file containing inital CDR region estimates in centromere. A tab delimited column (.tsv) withg chromosome as first column, start position of CDR region in second column and end position of CDR region of last column
- `-m` (default = False, optional) Determines method to use inout file for CDR regions initial transition and emission matrices or entered values. Default is using input file
- `-aa` (optional, float) Probability from next CpG position being in a CDR given current CpG position in CDR
- `-ab` (optional, float) Probability from next CpG position not being in a CDR given current CpG position in CDR
- `-ba` (optional, float) Probability from next CpG position not being in a CDR given current CpG position not in CDR
- `-bb` (optional, float) Probability from next CpG position being in a CDR given current CpG position not in CDR
- `-ax` (optional, float) Probability of current CpG position is methylated given current CpG position in CDR
- `-ay` (optional, float) Probability of current CpG position is not methylated given current CpG position in CDR
- `-bx` (optional, float) Probability of current CpG position is methylated given current CpG position not in CDR
- `-by` (optional, float) Probability of current CpG position is not methylated given current CpG position not in CDR


Description of method:

- Begins with an estimated path of methylation of CG sites (or any other methylated sites) within centromere regions
- Uses an EM algorithm using a first order HMM to determine CDRs
    -   E step: Estimate state path or if a each site is in a CDR or not in a CDR using Viterbi Algorithm
    -   M step: Resestimate transition and emission matrices using new state path and emission path aka Parameter Estimation
- Continues running until algorithm convergence/matrices stop changing
- Return bed file containing chromosome name, start position and end position of CDR estimates
