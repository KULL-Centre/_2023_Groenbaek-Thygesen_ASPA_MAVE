# _2023_Groenbaek-Thygesen_ASPA_MAVE
Scripts and output from "Deep mutational scanning reveals a tight correlation between protein degradation and toxicity of thousands of non-native aspartoacylase protein variants" by Martin Grønbæk-Thygesen, Vasileios Voutsinos, Kristoffer E. Johansson, Thea K. Schulze, Matteo Cagiada, Line Pedersen, Lene Clausen, Snehal Nariya, Rachel L. Powell, Amelie Stein, Douglas M. Fowler, Kresten Lindorff-Larsen, Rasmus Hartmann-Petersen.

Files and Directories

- abundance.csv: Abundance scores
- toxicity.csv: Toxicity scores
- tile_stability.csv: Degron TSI scores for ASPA tiles
- plots.ipynb: Jupyter notebook with plots
- pacbio: Scripts for processing long CCS reads into a barcode map 
- illumina: Scripts for calculating the abundance and toxicity scores

Pacbio and illumina sequencing reads for abundance and toxicity scores are available at https://doi.org/10.17894/ucph.3e05fe3a-4d7e-4d70-9056-18ed999e7e1e

Code used for degron scoring of tiles is available with the [Parkin VAMPseq paper by Clausen et al](https://github.com/KULL-Centre/_2023_Clausen_parkin_MAVE/tree/d9488dcfdb79329af2fae437d9f9452576d0a2d2/illumina_degron) with sequencing reads available at https://doi.org/10.17894/ucph.d879cfce-efb3-4eaa-928f-87a94d9560ef
