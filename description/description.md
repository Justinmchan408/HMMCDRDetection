# Method Description of HMM CDR

This markdown explains the methodology and the steps behind using the HMM Detection approach to identify chromosome dip regions, CDRs in centromere higher order repeat regions. CDR are regions of hypomethylation among regions of hypermethylation in highly repetetive alpha regions. These regions have been associated with sites of cenpA binding and kinetochore attachement.

## Disadvantages of Sliding Window Approach for CDR Detection

There are various methods that can be used to detect CDRs, one of them being the sliding window appraoch. The sliding window approach involves creating a window where all possible modification sites are observed. If a majority of sites are determined to be methylated, then this window is classified as a methylated region. If a majority of sites are determined to be not methylated, then this window is classified as a not methylated or CDR region. This may be helpful approach when exploring regions for the first time but, there are few disadvantages to the approach.

1. Different window sizes may be classify different regions differently

Some regions would be classified as both CDR Regions and methylated regions depending on the window sizes used in the sliding window algorithm. Since the size of interactions between kinetrochore proteins and centromere regions is unknown, it is difficult to estimate the ideal window size to be used when using the sliding window algorithm. Window sizes could possibly provide a probability of a maximum likelihood but, this would vary depending on the size of window sizes. A larger window size could give a more precise probability but, may eventually be too big to accuractly define insightful CDR regions while smaller windows would provide a less precise probability but, could provide where smaller CDR regions are located in centromere. 

2. Lack of Knowledge about Centromeres 

The lack of information on centromere regions and its proteins interactions is a difficult problem where a naive sliding window approach may not be the best method to accurately define CDR regions. Since little is known about the size of CDR regions, it would be a better approach to use a model to define CDRs rather than define a window length.

## Hidden Markov Model

A hidden Markov Model is a statistical Markov Model in which the system being modeled is assumed to be a Markov process and its states are hidden. Markov process is a stochastic process where there is a sequence of events where each event can be defined by a number of previous events in the sequence. For example, if an event was defined only by the previous event in a sequence, then it could be modeled by a first order Markov Model. If an event was defined by two previous events in a sequence, then it could be modeled by a second order Markov Models, etc. In a hidden Makov Model, the states are hidden where we do not know what phase or state of when each observation is occuring.

## Appling HMM to CDRs

A hidden Markov Model would be a great approach to detect CDRs since it fits the problem well. The problem in detecting CDRs in centromere involves two hidden states, CDR regions and non-CDR regions within the centromere. CDR regions would be locations in the centromere that proteins like cenpA are binding to while non-CDR regions would be methylated regions in the centromere. In addition to the states, there is observations which are possible modifications sites and their probabilities of a site being methylated. By using oservations of methylated and non-methylated sequence of sites, we can use the information to predict CDR regions

## Definining Transition and Emission Matrices

When implementing the Baum Welch Algorithm, it involves defining hidden states and emissions to run the algorithm. As mentioned previously, there can be two hidden states defined for the problem: CDR Regions and non-CDR regions that we cannot clearly observe from the data. There are two types of emissions for a possible modification site: methylated or non-methylated. Using the various different combinations of hidden states and emissions, we can create transition and emission matrices for the problem.

### Transition Matrix

The transition matrix is a matrix containing probabilities of the next event's hidden state given the current event's state. In other words, this is the probability that the next modification site is in a CDR or not in a CDR region, knowing that our current modification site is in a CDR or not in a CDR. The conditional probabilities can be written in the transition matrix that is used in the Baum-Welch Algorithm and other HMM related problems. The transition matrix is 2 x 2 matrix based on our 2 hidden states which can be defined as:

|                               | CDR (Next Mod Site)                |  Not a CDR (Next Mod Site)                |
| ----------------------------- | ---------------------------------- | ----------------------------------------- |
| CDR (Current Mod Site)        | -aa                                | -ab                                       |
| Not a CDR (Current Mod Site)  | -ba                                | -bb                                       |
