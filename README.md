# Epitope mapping of analytes binding RBD

This repository contains the tools used for the analysis of flow cytometry data to determine the epitopes of multiple analytes binding RBD.

## How to run the analysis
To install and run this project you need conda. Clone this repo to a location of your choice.
### Initial analysis of each residue position for each sample
1. Install the environment by running `conda env create -f environment.yml` from the main directory.
2. Acitvate the environment by running `conda activate capyflex-env`
3. Place the raw data in the data directory. Place all files belonging to a certain sample in a directory.
4. Run `code/capyflex_epitope.py -h` to see options.
5. From the main directory of your project run `code/capyflex_epitope.py -a <sample-name>` with the same sample name as the directory in 'data/'. Follow the instructons to run the script. If you are reproducing the analysis for the published project the `defaults.toml` files for each sample will be used.
6. Rerun step 5 for all samples.

### Plotting of epitopes
1. Run `code/capyflex_epiplot.py -h` to see options for the plotting script.
2. From the main directory of your project run `code/capyflex_epitope.py -p <positive-control> -n <negative-comntrol>` with the same control names as in your `.fcs` files. If you are reproducing the analysis for the problished project they will be `-p RBD -n FMO`.

### Finding the results
- Detailed flow cytometry plots for each sample in `results/fc_plots/`.
- Supplementary flow cytometry plots for each sample in `results/supplementary_fc_plots/`.
- Bar plots for each analyte and the combined heatmap for all analytes in `results/epitope_plots/`.
    - Values used to plot in `results/plot_values/`.

### NOTE
This repo is currently being set up and is not complete.