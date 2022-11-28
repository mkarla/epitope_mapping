#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:02:56 2022

@author: maximiliankarlander
"""

import argparse
from argparse import RawTextHelpFormatter
import tomlkit
import os
import glob
import sys
import math
from time import sleep
from tqdm import tqdm
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import cm
import seaborn as sns
import numpy as np
np.seterr(all="ignore")

"""
- pymol
    - pymol file located in pymol folder
    - show features
- inputparametrar
    - plot parameter
    - positive
    - negative
    - plot order
    - pymol on/off
    - pymol limits
    - structural
    - rotation
    - static areas + color
- defaults



- mappstruktur
- Heatmap
"""

def make_filelist(input_dir):
    if input_dir[-1] == '/':
        input_dir = input_dir
    else:
        input_dir = input_dir + '/' 
    filepath = input_dir  + 'output/results/'
    if not os.path.exists(filepath):
        print("The folder doesn't exist")
        sys.exit(1)
    files = glob.glob(filepath + '*.csv')     
        
    if len(files) == 0:
        print('No results files in input folder.')
        sys.exit(1)
    return(files)

def make_output(input_path):
    if input_path[-1] == '/':
        input_path = input_path
    else:
        input_path = input_path + '/' 
    output_path = input_path + 'output/'
    epitope_plot_path = output_path + 'epitope_plots/'
    if not os.path.exists(epitope_plot_path):
        os.makedirs(epitope_plot_path)
    pymol_path = output_path + 'pymol/'
    if not os.path.exists(pymol_path):
        os.makedirs(pymol_path)
    return(epitope_plot_path, pymol_path)

def NormalizeData(current):
    binding_control = current.loc[current["file"] == "RBD"]["value"].mean()
    no_binding_control = current.loc[current["file"] == "FMO"]["value"].mean()
    nomalized_data = [(x-no_binding_control)/(binding_control-no_binding_control)*100 for x in current["value"]]
    return(nomalized_data)

def data_collector(file, plot_parameter):
    df = pd.read_csv(file)
#    df = df[df["gated events"] > 100]
    df = pd.melt(df, id_vars=['file'], value_vars=df.columns[3:].values.tolist())
    df["file"] = [x.split("_")[0] for x in df["file"]]
    df = df.loc[df["variable"] == plot_parameter]
    df['normalized'] = NormalizeData(df)
    return(df)

def remove_replicates(df, positive):
    std_pos = df.loc[df["file"] == positive]["normalized"].values.tolist()
    std_pos = statistics.stdev(std_pos)
    df = df.groupby(['file']).mean()
    return(df, std_pos)

def order(df, order_param):
    if order_param == 'file':
        df = df.sort_values(by='file')
        plot_order = df.index.unique().tolist()
    elif order_param == 'value':
        df = df.sort_values(by='normalized', ascending=False)
        plot_order = df.index.unique().tolist()
    return(plot_order)

def plot_colors(df, value, vmin, vmax, structural):
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
    df = df[['log change']]
    df.columns = [analyte]
    df.index.names = ['Position']
    df2 = df.copy()
    df2.replace(-np.inf, 0, inplace=True)
    df2.replace(-np.nan, 0, inplace=True)
    return(df2)

def position_extractor(position):
    pos = ''
    for i in position:
        try: 
            k = int(i)
            if not type(k) == str:
                pos = pos+str(k)
        except:
            pass
    return(pos)
           
def pymol_scripter(path, summary, analyte, plot_parameter, increased, decreased, structural):
    no_change_color = summary.loc[summary["file"] == "RBD"]["rgba_colors"].values.tolist()[0]
    hex_col = matplotlib.colors.to_hex(no_change_color)
    hex_col = "0x" + hex_col[1:]
    
    #structural = summary.loc[summary['file'] == "RBD355"]
    filt_sum = summary.loc[(summary['file'] != "FMO") & (summary['file'] != "RBD355")]
    high_four = filt_sum.nlargest(increased, "change")
    low_four = filt_sum.nsmallest(decreased, "change")   
    extreme_sum = pd.concat([high_four, low_four, structural], axis=0)
    
    lines = [
    '#!/usr/bin/env python3',
    '# -*- coding: utf-8 -*',
    'from pymol import cmd',
    'cmd.delete("all")',
    'cmd.load("/users/maximiliankarlander/Dropbox/Jobb/Projects/epitope_mapping/pymol/6yla_annot.pse")',
    'cmd.bg_color(color="white")',
    'cmd.refresh()',
    f'cmd.color("{hex_col}", "rbd")',
    'cmd.color("0xFFFFFF", "rbd")',
    'cmd.hide("everything", "heavy")',
    'cmd.hide("everything", "light")',
    'cmd.rotate("y", "-110")',
    'cmd.color("0xa3a3a3", "resi 438-506")',
    'cmd.color("0x835A00", "resi 355")'
         ]

    if analyte == "CR3022":
        lines.append('cmd.show("cartoon", "L1")')
        lines.append('cmd.show("cartoon", "L2")')
        lines.append('cmd.show("cartoon", "L3")')
        lines.append('cmd.show("cartoon", "H1")')
        lines.append('cmd.show("cartoon", "H2")')
        lines.append('cmd.show("cartoon", "H3")')
        lines.append('cmd.color("0x8aac88", "cr3022_epitope_")')

    for pos, col in zip(extreme_sum["file"], extreme_sum["rgba_colors"]):  
            if position_extractor(pos) == 3:
                if pos in structural:
                    pos = position_extractor(pos)
#                    lines.append(f'cmd.color("0x549d47", "resi {pos}")') 
                    
                else:
                    pos = position_extractor(pos)
                    hex_col = matplotlib.colors.to_hex(col)
                    hex_col = "0x" + hex_col[1:]
                    lines.append(f'cmd.color("{hex_col}", "resi {pos}")')
    
    lines.append('cmd.refresh()')
    lines.append(f'cmd.png("{analyte}_structure_1.png", "30cm", "30cm", dpi=300, ray=1)')
    lines.append('cmd.turn("y", 180)')
    lines.append('cmd.refresh()')
    lines.append(f'cmd.png("{analyte}_structure_2.png", "30cm", "30cm", dpi=300, ray=1)')
        
                     
                     
    with open(path + analyte + "_6yla_" + plot_parameter + '_top' + str(increased) + '_bottom' + str(decreased) + '.py', 'w') as f:
        for l in lines:
            f.write(l + "\n")

def analysis_iterator(files, plot_parameter, barplot_path, pymol_path):
    positive = 'RBD'
    order_param = 'file'
#    plot_order = 'value'
    print('Determining the epitope for the following analytes:')
    files.sort()
    for file in files:
        print(file.split('/')[-1].split('.')[0])
    print('')
    print('Making individual plots for each analyte')
    for i in tqdm(range(len(files))):
        #ADD PYMOL
        file = files[i]
        sleep(0.1)
        analyte = file.split('/')[-1].split('.')[0]
        df = data_collector(files[i], plot_parameter)
        df = log_plot(df, plot_parameter, analyte, positive, order_param, barplot_path)
        df = heatmap_transform(df, analyte)
        if 'df_heatmap' not in locals():
            df_heatmap = df
        else: 
            df_heatmap = df_heatmap.join(df)
    df_heatmap = df_heatmap.T
    return(df_heatmap)

def plot_heatmap(df, heatmap_path):
    target = heatmap_path.split('/')[-4]
    cmap = sns.cm.crest_r
    cmap = 'Greens_r'
    fig, ax = plt.subplots(1,1, figsize=(20,5))
    ax = sns.heatmap(df, cmap=cmap, linewidth=.5, vmin=0.8, vmax=2.8)
    plt.tight_layout()
    df.to_csv('test.csv')
    plt.savefig(heatmap_path + target + '_epitopes.pdf', dpi=300)       

def main(args):
    input_dir = args.Input
    print('')
    print('Input directory:', input_dir)
    files = make_filelist(input_dir)
    
    plot_parameter = 'gated mean FL3-A div mean FL1-A one-by-one'
    epitope_plot_path, pymol_path = make_output(input_dir)
    
    print('')
    print('____________________________')
    print('Running')
    
    df_heatmap = analysis_iterator(files, plot_parameter, epitope_plot_path, pymol_path)
    
    print('')
    print('Making heatmap of all analytes')
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
#    parser.add_argument("-d", "--Defaults", help="Define new defaults", action='store_true')
    parser.add_argument("-i", "--Input", help="Input directory")
    parser.add_argument("-p", "--Positive", help="Name of positive control")
    parser.add_argument("-n", "--Negative", help="Name of negative control")
    parser.add_argument("-s", "--Struture", help="Name of negative control", action='store_true')
#    parser.add_argument("-hi", "--HighGate", help="Gate to remove non expressing population")
#    parser.add_argument("-dp", "--Double", help="!Experimental! Use this to gate out the lower population of two", action='store_true')
    args = parser.parse_args()  
    main(args)