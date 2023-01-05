# Epitope mapping of analytes binding RBD

This repository contains the tools used for the analysis of flow cytometry data to determine the epitopes of multiple analytes binding RBD.

## How to run the analysis
To install and run this project you need conda. Clone this repo to a location of your choice.
### Initial analysis of each residue position for each sample
The initial analysis is performed indepently for each analyte. This analysis reads all `.fcs` releated to one analyte and calculates the mean quotient of binding per expression for each sample (residue/control/etc). It will output relevant FC-plots, both detailed and more slimmed down, for each sample and a `.csv` file with all calculated values for that analyte.
1. Install the environment by running `conda env create -f environment.yml` from the main directory.
2. Acitvate the environment by running `conda activate capyflex-env`
3. Place the raw data in the data directory. Place all files belonging to a certain sample in a directory.
4. Run `code/capyflex_epitope.py -h` to see options.
5. From the main directory of your project run `code/capyflex_epitope.py -a <sample-name>` with the same sample name as the directory in 'data/'. Follow the instructons to run the script. If you are reproducing the analysis for the published project the `defaults.toml` files for each sample will be used.
6. Rerun step 5 for all samples.

### Plotting of epitopes
This part of the analysis takes the `.csv` files from the previous analysis steps, normalize the data, and make both individual plots for each analyte and a combined heatmap. Furhtermore there are options to make a PyMol script to show the top three disrupting positions in a pdb structure.
1. Run `code/capyflex_epiplot.py -h` to see options for the plotting script.
2. From the main directory of your project run `code/capyflex_epiplot.py -p <positive-control> -n <negative-comntrol>` with the same control names as in your `.fcs` files. If you are reproducing the analysis for the published project they will be `-p RBD -n FMO`.
3. To write PyMol scripts add the `-py`command. This will furhther require the `-f <framework-name>` command where the name of the framework in the PyMol file used. The PyMol can either be `.pdb` or `.pse` with the later having priority. Further options are `-r <initial-rotation>`, `-s <structural-positions>`, and `-hi <regions-to-hide>`. When running the PyMol script it will color the framework region white, structural positions brown, and the top three disrupting positions colored red with intensity correlating to its binding value. FInally two `.png` files be saved, the second being rotated 180° on the y axis.

### Finding the results
- Detailed flow cytometry plots for each sample in `results/fc_plots/`.
- Supplementary flow cytometry plots for each sample in `results/supplementary_fc_plots/`.
- Bar plots for each analyte and the combined heatmap for all analytes in `results/epitope_plots/`.
    - Values used to plot in `results/plot_values/`.
- PyMol scripts in `results/pymol`.

### NOTE
This repo is currently being set up and is not complete.

### Stuff of importance
Following plot command:
`code/capyflex_epiplot.py -n FMO -p RBD -py -f rbd -r -110 -s RBD355R -hi heavy light`