How to convert sci-Hi-C \*.matrix files to cell-LP matrix (./data_conversion):
1. Run matrix_to_interaction.py to convert the sci-Hi-C data to FitHiC format
2. Run interaction_to_cisTopic_sparse.py to convert to sparse matrix format (cell_idx | LP_idx | count)
3. Run concatenate_sparse_mat.sh to combine sparse contact matrices into one file
4. Run concatenate_samples_sparse_mat.sh to combine concatenated files from Step 4 into one big file
        Run this if you have multiple samples/libraries
5. Run run_cisTopic_sparse.R to perform topic modeling


