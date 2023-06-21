# O-GlcNAcPRED-DL-master

## Introduction
O-GlcNAcPRED-DL is a tool to predict protein O-GlcNAcylation sites on an ensemble model of deep learning.

Webserver is online, which is freely accessible at https://oglcnac.org/pred_dl/


## Directory
[feature] Used for encoding, which contain AAindex encoding method and Word2Vec encoding method. 

[model] Models trained for Human and Mouse respectively(namely O-GlcNAcPRED-DL). 

[seq_cut] Given a protein sequence, peptide segments with a growth degree of 29 can be symmetrically truncated left and right centered around S/T.

[test_data] The independent test set used in this experiment and the data of Arabidopsis, Caenorhabditis elegans, Drosophila, Rat, Wheat, and Others species in the cross species experiment.


## Usage
Hwmean_voting.py: predict Human protein O-GlcNAcylation sites.

Mwmean_voting.py: predict Mouse protein O-GlcNAcylation sites.
