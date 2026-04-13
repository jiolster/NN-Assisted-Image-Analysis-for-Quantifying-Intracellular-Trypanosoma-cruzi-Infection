# NN-Assisted Image Analysis for Quantifying Intracellular Trypanosoma cruzi Infection
Repository of the scripts used for the publication "NN-Assisted Image Analysis for Quantifying Intracellular Trypanosoma cruzi Infection", at the Signaling and Adaptive Mechanisms in Trypanosomatids Laboratory.

The main scripts can be run on Google Colab with ContandoAmas.ipynb (Español) or CountingAmastigotes.ipynb (English).

The models that were trained for this work are available on Hugging Face:
https://huggingface.co/jiolster/Automated-Image-Analysis-of-the-Intracellular-Stage-of-Trypanosoma-cruzi/tree/main

Note: the main data used for the manuscript were analyzed using the RecuentoNN.py script, which generated the sc_results.csv and fov_results.csv tables (inside Statistics and figures/Data). This script was written to use images contained in seperate folders that contain all channels for each field of view. 

Other data used for the manuscritp: data containing "manual correspond to the manual quantification, data with "inspect" refers to quantifications obtained using the iamge analyisis methodology replciated from Yazdanparast et al. (2014).
