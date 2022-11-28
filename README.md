# Epitope mapping of analytes binding RBD

This repository contains the tools used for the analysis of flow cytometry data to determine the epitopes of multiple analytes binding RBD.

## How to run the analysis
To install and run this project you need conda. Clone this repo to a location of your choice.
1. Install the environment by running `conda env create -f environment.yml` from the main directory.
2. Acitvate the environment by running `conda activate capyflex-env`
3. cd to the cloned directory and place the raw data in the data directory. Place all files belonging to a certain sample in a directory.
4. Run `python3 code/capyflex -h` to see options.
4. From the main directory of your project run `python3 code/capyflex_epitope -a <sample-name>` with the same sample name as the directory in 'data/'. Follow the instructons to run the script.

### NOTE
This repo is currently being set up and is not complete.