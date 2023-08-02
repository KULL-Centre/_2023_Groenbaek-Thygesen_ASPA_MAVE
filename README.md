# _2023_Groenbaek-Thygesen_ASPA_MAVE
Scripts and output from "Deep mutational scanning reveals a tight correlation between protein degradation and toxicity of thousands of non-native aspartoacylase protein variants" by Martin Grønbæk-Thygesen, Vasileios Voutsinos, Kristoffer E. Johansson, Thea K. Schulze, Matteo Cagiada, Line Pedersen, Lene Clausen, Snehal Nariya, Rachel L. Powell, Amelie Stein, Douglas M. Fowler, Kresten Lindorff-Larsen, Rasmus Hartmann-Petersen.

Files and Directories

- abundance.csv: Abundance scores
- toxicity.csv: Toxicity scores
- tile_stability.csv: Degron TSI scores for ASPA tiles
- plots.ipynb: Jupyter notebook with plots based on the article supplementary CSV file
- pacbio: Scripts for processing long CCS reads into a barcode map 
- illumina: Scripts for calculating the abundance scores
- illumina_degron: Scripts for calculating the degron TSI scores

Code used for degron scoring of tiles is available with the [Parkin VAMPseq paper by Clausen et al](https://github.com/KULL-Centre/_2023_Clausen_parkin_MAVE/tree/d9488dcfdb79329af2fae437d9f9452576d0a2d2/illumina_degron)

Pacbio and illumina sequencing reads for abundance scores are available at 

Sequencing reads for degron TSI scores are available at https://doi.org/10.17894/ucph.d879cfce-efb3-4eaa-928f-87a94d9560ef
