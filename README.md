# Capturing cell-type specific compartment patterns by applying topic modeling to single-cell Hi-C data

This is repository for the manuscript "Capturing cell-type specific compartment patterns by applying topic modeling to single-cell Hi-C data". 

# Data

Raw .fastq and aligned .bam files are available on 4DN Data Portal (https://data.4dnucleome.org/). <br />
Processed data files are available on https://noble.gs.washington.edu/proj/schic-topic-model.

# Data preprocessing
For preprocessing our data we followed the pipeline in https://github.com/VRam142/combinatorialHiC to align the reads to hg19 and generate sci-Hi-C matrix files binned at 500kb resolution. The .matrix files are in the format:
> bin1  bin2  count normalized_count  chr1  chr2

# Construction of cell-locus pair matrix for topic modeling
1. Run matrix_to_interaction.py to convert sci-Hi-C .matrix files to FitHiC interaction format (.int.bed):
> chr1  midpoint1 chr2  midpoint2 count normalized_count
2. Run interaction_to_sparse_matrix.py to convert FitHiC interaction files to sparse matrix format: 
> cell_idx | LP_idx | count
3. Run concatenate_sparse_mat.sh to combine sparse contact matrices into one file
4. Run concatenate_samples_sparse_mat.sh to combine concatenated files from Step 4 into one big cell-LP matrix. This code is useful for combining matrices from multiple samples/libraries.
5. Run run_cisTopic_sparse.R to perform topic modeling on the cell-locus pair matrix.
