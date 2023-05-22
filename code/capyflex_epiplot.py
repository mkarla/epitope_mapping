#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:02:56 2022

@author: maximiliankarlander
"""

import argparse
from argparse import RawTextHelpFormatter
import os
import glob
import sys
import math
from tqdm import tqdm
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import cm
import seaborn as sns
import numpy as np
np.seterr(all="ignore")

def make_filelist(input_dir, args):
    #This function makes a list of all sample summaries to include
    if input_dir[-1] == '/':
        input_dir = input_dir
    else:
        input_dir = input_dir + '/' 
    filepath = input_dir
    if not os.path.exists(filepath):
        print("The folder doesn't exist")
        sys.exit(1)
    files = glob.glob(filepath + '*.csv')     
        
    if len(files) == 0:
        print('No results files in input folder.')
        sys.exit(1)
    
    if args.Custom:
        files, tops = custom_files(args.Custom, files)
    else:
        tops=[args.Top]*len(files)
    return(files, tops)

def pymol_check(args):
    #This function checks if pymol script should be written and if the correct parameters have been given
    if args.Pymol:
        print('Also making PyMol scripts for each analyte')
        if not args.Framework:
            print('      - No framework provided. Please provide the name of the framework (-f) in model.')
            sys.exit(1)
        if args.Structural:
            print('      - The following residues will be considered structurally important and not considered for the epitope:')
            for residue in args.Structural:
                print(f'         - {residue}')

def path_check(path):
    if not os.path.exists(path):
        os.makedirs(path)

def make_output():
    #This function makes output directories
    output_path = 'results/'
    epitope_plot_path = output_path + 'epitope_plots/'
    path_check(epitope_plot_path)
    pymol_path = output_path + 'pymol/'
    path_check(pymol_path)
    plot_value_path = output_path + 'plot_values/'
    path_check(plot_value_path)
    return(epitope_plot_path, pymol_path, plot_value_path)

def NormalizeData(current, positive, negative):
    #This function normalizes the data
    binding_control = current.loc[current["file"] == positive]["value"].mean()
    no_binding_control = current.loc[current["file"] == negative]["value"].mean()
    nomalized_data = [(x-no_binding_control)/(binding_control-no_binding_control)*100 for x in current["value"]]
    return(nomalized_data)

def data_collector(file, plot_parameter, positive, negative):
    #This function reads the data of a summary, reformats it, and normalize it
    df = pd.read_csv(file)
    df = pd.melt(df, id_vars=['file'], value_vars=df.columns[3:].values.tolist())
    df["file"] = [x.split("_")[0] for x in df["file"]]
    df = df.loc[df["variable"] == plot_parameter]
    df['normalized'] = NormalizeData(df, positive, negative)
    return(df)

def remove_replicates(df, positive):
    #This function replaces replicates with the mean values and provides the std
    std_pos = df.loc[df["file"] == positive]["normalized"].values.tolist()
    std_pos = statistics.stdev(std_pos)
    df = df.groupby(['file']).mean('value', 'normalized')
    return(df, std_pos)

def order(df, order_param):
    #This function order plot data based on how the user want to order it, either by file name or value
    if order_param == 'file':
        df = df.sort_values(by='file')
        plot_order = df.index.unique().tolist()
    elif order_param == 'value':
        df = df.sort_values(by='normalized', ascending=False)
        plot_order = df.index.unique().tolist()
    return(plot_order)

def plot_colors(df, value, vmin, vmax, structural):
    #This function assigns a color to each value for plotting based on the value
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)    
    rgba_colors = []
    positions = df.index.unique().tolist()
    for i in positions:
        position = df.loc[df.index == i]
        if i in structural:
            rgba_color = (0.5137254902, 0.3529411765, 0, 1)
        else:
            if math.isnan(position.iloc[0][value]):
                rgba_color = (0.705673158, 0.01555616, 0.150232812, 1.0)
            else:
                rgba_color = cm.coolwarm_r(norm(position.iloc[0][value]))
        rgba_colors.append(rgba_color) 
    df["rgba_colors"] = rgba_colors  
    return(df, rgba_colors)

def log_plot(df, plot_paramater, analyte, positive, order_param, barplot_path):
    #This function makes a bar plot of all values for an analyte
    df, std_pos = remove_replicates(df, positive)
    df["log change"]=np.log10(df["normalized"])
    plot_order = order(df, order_param)
    df, rgba_colors = plot_colors(df, 'log change', 0.8, 3.2, [])
    fig, ax = plt.subplots(1,1, figsize=(20,5)) 
    ax = sns.barplot(x=df.index, y="normalized", data=df, order=plot_order, palette=rgba_colors, ax=ax, log=True)
    ax.set_ylabel("normalized binding")
    ax.set_ylim([1,1000]) 
    ax.set_xlabel("Alanine mutant")
    ax.set_xlim(-1, len(df))
    ax.set_title(analyte, size=20)
    plt.xticks(rotation = 90)
    plt.tight_layout()
    plt.savefig(barplot_path + analyte + '.pdf', dpi=300)
    return(df)

def heatmap_transform(df, analyte):
    #This function transform the data for each analyte for making the heatmap
    df = df[['log change']]
    df.columns = [analyte]
    df.index.names = ['Position']
    df2 = df.copy()
    df2.replace(-np.inf, 0, inplace=True)
    df2.replace(-np.nan, 0, inplace=True)
    return(df2)

def position_extractor(position):
    #This function extracts the mutated position from the filename
    pos = ''
    for i in position:
        try: 
            k = int(i)
            if not type(k) == str:
                pos = pos+str(k)
        except:
            pass
    return(pos)

def custom_files(custom, files):
    #This function handles custom input of what analytes to be analyzed and how many top residues to show in pymol scripts
    column = 0
    custom_list = []
    for i in custom:
        if column==0:
            current = [i]
            column=1
        elif column==1:
            current.append(i)
            custom_list.append(current)
            column= 0
    custom_list.sort(key=lambda x:x[0])
    custom_files = []
    tops = []
    for i in custom_list:
        custom_file = False
        top = False
        for file in files:
            if i[0] == file.split('/')[-1].split('.')[0]:
                custom_file = file
                try:
                    top = int(i[1])
                except:
                    print("Non correct input of -c. Input should be each analyte followed by desired number of top residues. For example -c MO176_156 5 MO176_301 2")
                    sys.exit(1)
        if custom_file and top and len(custom)%2==0:
            custom_files.append(custom_file)
            tops.append(top)
        else:
            print("Non correct input of -c. Input should be each analyte followed by desired number of top residues. For example -c MO176_156 5 MO176_301 2")
            sys.exit(1)
    return(custom_files, tops)

def pymol_scripter(analyte, path, df, neg, struct, rotation, framework, hide, top):
    #locate pdb and pse files. pse file have priority over pdb.
    pdb_file = glob.glob("data/pymol/" + '*.pdb')
    pse_file = glob.glob("data/pymol/" + '*.pse')
    if len(pse_file)==1:
        file = os.path.abspath(pse_file[0])
    elif len(pdb_file)==1:
        file = os.path.abspath(pdb_file[0])
    else:
        print('Too many pdb or pse files were provided.')
        sys.exit(1)
    model = file.split("/")[-1]

    #initiate lines for pymol script
    lines = ['#!/usr/bin/env python3',
            '# -*- coding: utf-8 -*',
            '',
            '#Initialize PyMol and load file',
            'from pymol import cmd',
            'cmd.delete("all")',
            f'cmd.load("{file}")',
            'cmd.bg_color(color="white")',
            'cmd.refresh()']

    #define frameworks and regions to hide
    lines.append('')
    lines.append("#Regions to show and hide, including structural residues")
    for region in framework:
        lines.append(f'cmd.color("0xFFFFFF", "{region}")')
    if hide:
        for region in hide:
            lines.append(f'cmd.hide("everything", "{region}")')
    if struct:
        for residue in struct:
            df = df.loc[df.index!=residue]
            residue = position_extractor(residue)
            lines.append(f'cmd.color("0x835A00", "resi {residue}")')
    
    #color top residues according to binding change
    lines.append('')
    lines.append("#Top residues")
    
    df = df.loc[df.index!=neg]
    """
    df_top =df.nsmallest(top, "log change")
    for pos, col in zip(df_top.index, df_top["rgba_colors"]):   
        pos = position_extractor(pos)
        hex_col = matplotlib.colors.to_hex(col)
        hex_col = "0x" + hex_col[1:]
        lines.append(f'cmd.color("{hex_col}", "resi {pos}")') """  
    df = df.sort_values(by='log change')
    count = 0
    for pos, col in zip(df.index, df["rgba_colors"]):
        pos = position_extractor(pos)
        hex_col = matplotlib.colors.to_hex(col)
        hex_col = "0x" + hex_col[1:]
        if count < top:
            lines.append(f'cmd.color("{hex_col}", "resi {pos}")')
        count += 1

    #rotation and save
    lines.append('')
    lines.append("#Refresh and save .png files from two angles")
    if rotation:
        rotation_init = f'{rotation}'
        lines.append(f'cmd.rotate("y", "{rotation_init}")')
    lines.append('cmd.refresh()')
    lines.append(f'cmd.png("{analyte}_top={top}_structure_1.png", "30cm", "30cm", dpi=300, ray=1)')
    lines.append('cmd.turn("y", 180)')
    lines.append('cmd.refresh()')
    lines.append(f'cmd.png("{analyte}_top={top}_structure_2.png", "30cm", "30cm", dpi=300, ray=1)')

    #All residues
    lines.append('')
    lines.append("#All residues with corresponding color")
    for pos, col in zip(df.index, df["rgba_colors"]):
        pos = position_extractor(pos)
        hex_col = matplotlib.colors.to_hex(col)
        hex_col = "0x" + hex_col[1:]
        if pos != "":
            lines.append(f'#cmd.color("{hex_col}", "resi {pos}")')

    with open(f'{path}{analyte}_{model}_top={top}.py', 'w') as f:
        for l in lines:
            f.write(l + "\n")

def analysis_iterator(files, plot_parameter, barplot_path, pymol_path, plot_value_path, args, tops):
    #This function itereates over each analyte summary and collects processed data
    positive = args.Positive
    negative = args.Negative
    order_param = 'file'
    print('Determining the epitope for the following analytes:')
    files.sort()
    tab_width = 7
    for file in files:
        name = file.split('/')[-1].split('.')[0]
        if len(name) > tab_width:
            tab_width = len(name)
    tab_width = tab_width + 3
    print('Analyte', " "*(tab_width-7), "No of residues for PyMol (if used)")
    for i in range(len(files)):
        name = files[i].split('/')[-1].split('.')[0]
        print(name, " "*(tab_width-len(name)), tops[i])
    print('\nMaking individual plots for each analyte')
    pymol_check(args)
    for i in tqdm(range(len(files))):
        file = files[i]
        analyte = file.split('/')[-1].split('.')[0]
        df = data_collector(files[i], plot_parameter, positive, negative)
        df.to_csv(plot_value_path + analyte + '.csv')
        df = log_plot(df, plot_parameter, analyte, positive, order_param, barplot_path)
        if args.Pymol:
            pymol_scripter(analyte, pymol_path, df, negative, args.Structural, args.Rotation, args.Framework, args.Hide, tops[i])
        df = heatmap_transform(df, analyte)
        if 'df_heatmap' not in locals():
            df_heatmap = df
        else: 
            df_heatmap = df_heatmap.join(df)
    df_heatmap = df_heatmap.T
    return(df_heatmap)

def plot_heatmap(df, heatmap_path):
    #This function plots the heatmap
    height = 2.5 + len(df)/2
    name = list(df.index.unique())
    name = ('_').join(name)
    cmap = sns.cm.crest_r
    cmap = cm.coolwarm_r
    fig, ax = plt.subplots(1,1, figsize=(30,height))
    ax = sns.heatmap(df, cmap=cmap, linewidth=.5, vmin=0.8, vmax=3.2)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', rotation=0, labelsize=24)
    plt.xlabel('Position', size=24)
    plt.tight_layout()
    plt.savefig(f'{heatmap_path}epitope_heatmap_{name}.pdf', dpi=300)       

def main(args):
    input_dir = args.Input
    print(f'\nInput directory: {input_dir}') 
    files, tops = make_filelist(input_dir, args)
    plot_parameter = 'gated mean FL3-A div mean FL1-A one-by-one'
    epitope_plot_path, pymol_path, plot_value_path = make_output()
    
    print('\n____________________________\nRunning')
    
    df_heatmap = analysis_iterator(files, plot_parameter, epitope_plot_path, pymol_path, plot_value_path, args, tops)
    print('\nMaking heatmap of all analytes')
    plot_heatmap(df_heatmap, epitope_plot_path)
       
if __name__ == '__main__':
    prog = "Epitope mapping plotting program"
    description = """Version 0.1, author Maximilian Karlander \n
    This script is designed to take the output files of multiple analytes from capyflex_epitope and reveal their epitopes.
    """
    epilog = "I did good, didn't I?"
    parser = argparse.ArgumentParser(prog=prog, 
                                     description=description, 
                                     epilog=epilog,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--Input", help="Input directory", nargs='?', const="results/csv/", type=str, default="results/csv/")
    parser.add_argument("-p", "--Positive", help="Name of positive control")
    parser.add_argument("-n", "--Negative", help="Name of negative control")
    parser.add_argument("-py", "--Pymol", help="Make pymol script", action='store_true')
    parser.add_argument("-f", "--Framework", nargs='+', help="Name of framework region(s) for pymol script")
    parser.add_argument("-hi", "--Hide", nargs='+', help="Regions to hide for pymol script")
    parser.add_argument("-r", "--Rotation", help="Pymol initial rotation on y axis")
    parser.add_argument("-s", "--Structural", nargs='+', help="Structural positions")
    parser.add_argument("-t", "--Top", help="Number of positions to highlight in pymol script. Default is 3", nargs='?', const=3, type=int, default=3)
    parser.add_argument("-c", "--Custom", nargs='+', help="Analytes to include and number of positions to highlight in pymol script per analyte. This feature overrides --Top.")
    args = parser.parse_args()  
    main(args)