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
This part of the analysis takes the `.csv` files from the previous analysis steps, normalize the data, and make both individual plots for each analyte and a combined heatmap. Furhtermore there are options to make a PyMOL script to show the top disrupting positions in a pdb structure.
1. Run `code/capyflex_epiplot.py -h` to see options for the plotting script.
2. From the main directory of your project run `code/capyflex_epiplot.py -p <positive-control> -n <negative-comntrol>` with the same control names as in your `.fcs` files. If you are reproducing the analysis for the published project they will be `-p RBD -n FMO`.
3. To write PyMOL scripts add the `-py`command. Se below for detailed instructions.

#### Making PyMOL scripts
Part of the plotting includes the possibility to make a PyMOL script which will highlight a certain number of residues in a pdb file with a color correspoding to its heatmap color. The script will then save two images with 180° rotation of your protein. There are several options for doing this which will be outlined below.
- Your pdb file can be either a raw `.pdb` file, or you can use an annotated and adjusted `.pse` file. While the software is agnostic of your choice the latter is recommended as it will grant you more control of the output.
- Assuming you're using a `.pse` file, make a named selection of your framework in PyMOL. This will later be used by the script. Make sure only this region and regions you may want to use later (as named selections) are visible.
- To make a PyMOL script add the `-py` tag when running the plotting script. It will require the following tag as well:
    - `-f <framework-name>` command where the name of the framework in the PyMOL file used.
    - This will generate a PyMOL script where the top three positions are shown
- Apart from the required tags there are several optional tags:
    - `-hi <regions-to-hide>` This option lets you hide named selections you want to hide.
    - `-r <initial-rotation>` This option lets you make an initial rotation on the y axis.
    - `-s <structural-positions>` This option lets you specify residues you don't want to include in your determination of the epitope. One such reason may be if you identify a residue which may be of structural importance for the protein.
    - `-t <number-of-top-residues>` This option lets you choose how many residue to show in the images. The default option is `3`.
    - `-c <subset-of-analytes>` This option lets you choose a subset of you analytes to look at. Your chosen subset will also be reflected in the heatmap. Apart from only specifiying your analytes you will also need to specify how many residues to show for each analyte. This gives you the option to choose a different number of residues for different analytes. 
        - Example: you have three analytes called `one`, `two`, and `three` but only want to include `one` and `three`. For `one` you want to show the top three residues and for `three` you want to show the top five residues. In this case you would add `-c one 3 three 5`.
        - Note that using the `-c` tag overrides the `-t` tag.

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