# HFOUDCA

This folder contains:
Analysis of EEG data from EDF file format. Inlcudes 60 line noise removal,
surface laplacian, 3 (and growing) different detectors, and unsupervised
diffusion component analysis.

--------------------------------------------------------------------------
Parameters:
Out:
V - Eigenvectors with 20 largest eigenvalues (listed in descending order)
D - Respective 20 largest eigenvalues
HFO_Count - HFO_Count(i) is the number of HFO detections in channel i.
Therefore, sum(HFO_Count)=size(V,1)
filtData - Filtered detected HFOs with baseline removed.
rawData - Raw detected HFO data.
set - Lists start and stop times of all detections
info - Detector info to help identify files.
In:
path - File path of EDF file to be analyzed
locPath - File path of channels to use. Equivalent of EEG.chanlocs
eegPath - Path to folder that contains eeglab
chan - index of channels to use ex:[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 19 29 30]
detect - detector to be used. Currently allows "STE", "HIL", and "MNI"

--------------------------------------------------------------------------
NOTE: Currently, surface Laplacian function is from MikeXCohen at
http://mikexcohen.com and detectors are modified from RippleLab.
Otherwise, entirely coded by Jiaju Liu at jiajuliu1@gmail.com as part of
EPIBIOS4Rx project.
