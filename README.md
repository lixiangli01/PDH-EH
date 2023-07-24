# PDH-EH

Intro
-
PDH-EH is a novel framework for predicting DNA-protein binding hotspots using fused features. It leverages both embedded features extracted from protein language pre-trained models and handcrafted features for prediction.

Usage
-
We offer two methods for predicting hotspot residues. One is using fused features for prediction, while the other uses only embedded features for prediction.

1、Feature extraction.

For prediction based on fused features, it is necessary to extract manual features using the software below as input for the program. The extraction method is as follows:

Amino acid physicochemical properties: Retrieved from aaindex.txt.

Solvent accessibility-related features: Obtained from the NACCESS program.

Electrostatic potential-related features: Obtained from the APBS program.

Hydrogen bond-related features: Obtained from the HBPLUS program.

Secondary structure-related features: Obtained from the DSSP program.

Position-specific scoring matrix (PSSM) related features: Obtained from the PSI-BLAST program.

spatial neighbor-based PSSM (SNB-PSSM) related features: Obtained from the SNB-PSSM program.

Please save the above features (a total of 117 dimensions, while the remaining embedded features will be generated automatically in the prediction program) in a CSV format file, for example: 1CKTS41.csv. Additionally, please download the Bert model ("pytorch_model.bin") from [ProtTrans](https://huggingface.co/Rostlab) and place it in the "set_file" folder.

2、Prediction.

If you want to predict the amino acid S41 in the 1CKT complex, use the following command for prediction:
```
python predict_EH.py 1CKT A S 41 1CKTS41.csv
```

If you want to predict using only embedded features, use the following command:
```
python predict_E.py 1CKT A S 41 
```
The program's output result is the probability of predicting the residue as a hotspot.
