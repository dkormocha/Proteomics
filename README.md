# Proteomics
MSc Bioinformatics Project: Explaining incomplete tryptic digestion

The aim of the project was to explain incomplete tryptic digestion and create machine learning models to predict missed cleavage sites.
This repsoritory contains all the code used for data wrangling of mass spectrometry proteomics data plus, building Random Forest and Neural Network models. 


MSc Bioinformatics| Dhiraj Kormocha

# Explaining incomplete tryptic digestion
By Dhiraj Kormocha
Supervisor: Professor Conrad Bessant
2022

This research dissertation is submitted for the MSc in Bioinformatics at Queen Mary, University of London
MSc Bioinformatics| Dhiraj Kormocha

## Abstract
Trypsin is a highly specific proteolytic enzyme that exclusively hydrolyses peptide bonds at the carboxyl-terminus of Lysine and Arginine. This specificity makes trypsin the preferred proteolytic enzyme for protein digestion in bottom-up proteomics. This process requires trypsin to completely digest proteins in a complex sample into small peptides to identify and quantify proteins with high confidence. Unfortunately, during this process proteins are not completely digested, this is known as incomplete tryptic digestion and as a consequence missed cleaved peptides are produced. Predicting these missed cleavage sites can be a challenging problem, however, in this research machine learning models such as random forest and neural networks were used to predict missed cleavage sites. The existing missed cleavage predictor only uses the amino acid sequence and does not consider the 3D structure of the protein. This research mainly focuses on the 3D structure and its properties to predict missed cleavage. Features such as secondary structure and residue depth were incorporated into the machine learning models. The outcome of this project was a random forest model with a 0.84 accuracy and 0.84 precision in predicting missed cleavage; although, the neural network model performed poorly with a low accuracy possibly due to weak optimisation and small training dataset.

Acknowledgments:
A special thank you to Professor Conrad Bessant for the guidance and support during this project and throughout the bioinformatics course. Thank you to Esteban Gea and the Bessant lab for the support and resources provided during the project.
