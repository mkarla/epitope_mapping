#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 09:26:38 2022

@author: maximiliankarlander
"""
import argparse
from argparse import RawTextHelpFormatter
import tomlkit
import os
import glob
import sys
from scipy.optimize import curve_fit
import matplotlib.patches as mpl_patches
import FlowCal
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import statistics
from time import sleep
from tqdm import tqdm
np.seterr(all="ignore")

"""
- dela upp för snakemake
- generalisera adaptive ellipse?
"""

def argument_check(args):
    if not args.Input:
        print('No input path given, use -i to define input path')
        sys.exit(1)
    if not args.Analyte:
        print('No analyte given, use -a to define analyte')
        sys.exit(1)

def make_defaults(args, files, input_path, analyte):
    print('')
    savename = input_path + analyte + '/' + 'defaults.toml'
    toml = tomlkit.document()
    toml.add(tomlkit.comment('Defaults file for epitope mapping'))
    toml.add(tomlkit.nl())    
    center, channels = center_find(files, input_path)
    toml.add('center_x', center[0])
    toml.add('center_y', center[1])
    if not args.Channel1:
        print('No fluorescent channel given for expression, use -ch1 to define the channel')
        print('Available channels:', channels)
        channel1 = input('Channel: ')
    else:
        channel1 = args.Channel1
    while channel1 not in channels:
        print('Invalid channel 1.')
        channel1 = input('Channel: ')
    if not args.Channel2:
        print('No fluorescent channel given for binding, use -ch2 to define the channel')
        print('Available channels:', channels)
        channel2 = input('Channel: ')
    else:
        channel2 = args.Channel2
    while channel2 not in channels:
        print('Invalid channel 2:')
        channel2 = input('Channel: ')
    if not args.HighGate:
        print('No high gate given for expression, use -hi to define the gate')
        high_gate = input('Cutoff: ')
    else:
        high_gate = args.HighGate
    adapt_ellipse = args.Double 
    
    toml.add('channel1', channel1)
    toml.add('channel2', channel2)
    toml.add('high_gate', high_gate)
    toml.add('adapt_ellipse', adapt_ellipse)
    toml.add('channels', channels)
    with open(savename, mode="wt", encoding="utf-8") as fp:
        tomlkit.dump(toml, fp)
    return(center, channel1, channel2, high_gate, adapt_ellipse, channels)

def parameters(args, input_path, analyte, files):
    if input_path[-1] == '/':
        input_path = input_path
    else:
        input_path = input_path + '/'
    defaults_name = input_path + analyte + '/' + 'defaults.toml'
    if args.Defaults:
        print('Overwriting defaults')
        center, channel1, channel2, high_gate, adapt_ellipse, channels = make_defaults(args, files, input_path, analyte)
    elif os.path.exists(defaults_name):
        print('Using previous settings from defaults.')
        defaults_name = input_path + analyte + '/' + 'defaults.toml'
        with open(defaults_name, mode="rt", encoding="utf-8") as fp:
            defaults = tomlkit.load(fp)
        center = (defaults['center_x'], defaults['center_y'])
        channel1 = defaults['channel1']
        channel2 = defaults['channel2']
        channels = defaults['channels']   
        high_gate = defaults['high_gate']
        adapt_ellipse = defaults['adapt_ellipse']
    else:
        print('No defaults file found, manual input needed. Your input will be saved to a default file.')
        print('To overwrite the default file in the future, use the -d command')
        center, channel1, channel2, high_gate, adapt_ellipse, channels = make_defaults(args, files, input_path, analyte)
    if args.Channel1:
        channel1 = args.Channel1
    if args.Channel2:
        channel2 = args.Channel2
    while channel1 not in channels:
        print('Invalid channel 1.')
        channel1 = input('Channel: ')
    while channel2 not in channels:
        print('Invalid channel 2.')
        channel2 = input('Channel: ') 
    if args.HighGate:
        high_gate = args.HighGate
    if args.Double:
        adapt_ellipse = True
    return(center, channel1, channel2, high_gate, adapt_ellipse)

def make_filelist(input_dir, analyte):
    if input_dir[-1] == '/':
        input_dir = input_dir
    else:
        input_dir = input_dir + '/' 
    filepath = input_dir + analyte + '/'
    if not os.path.exists(filepath):
        print("The folder doesn't exist")
        sys.exit(1)
    files = glob.glob(filepath + '*.fcs')     
        
    if len(files) == 0:
        print('No FCS files in input folder.')
        sys.exit(1)
    return(files)

def make_output(analyte):
    output_path = 'results/'
    fc_plot_path = output_path + 'fc_plots/' + analyte + '/'
    if not os.path.exists(fc_plot_path):
        os.makedirs(fc_plot_path)
    sup_fc_plot_path = output_path + 'supplementary_fc_plots/' + analyte + '/'
    if not os.path.exists(sup_fc_plot_path):
        os.makedirs(sup_fc_plot_path)
    results_path = output_path + 'csv/'
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    return(fc_plot_path, sup_fc_plot_path, results_path)

def results_format(channel1, channel2):
    results = {
        "file": [],
        "gated events": [],
        "gated mean "+channel2: [],
        "gated mean "+channel1: [],
        "gated mean "+channel2+" div mean "+channel1: [],
        "gated mean "+channel2+" div mean "+channel1+" one-by-one": [],
        "gated mean "+channel2+" div mean "+channel1+" one-by-one stdev": []
        }
    return(results)

def center_find(filelist, output_path):
    print('Centertest being performed on:', filelist[0])
    s, s_g1, events, gates = read_enhance_fcs(filelist[0])
    df = pd.DataFrame(data=s_g1, columns=s.channels)
    meanSSC = df['SSC-A'].mean()
    meanFSC = df['FSC-A'].mean()
    center = (math.log(meanSSC, 10), math.log(meanFSC, 10))
    print('Suggested center (mean SSC, mean FSC):', center)
    done = 'n'
    channels = s.channels
    while done != 'y':
        s, s_g1, events, gates = read_enhance_fcs(filelist[0])
        fig, ((ax1), (tb1)) = plt.subplots(2,1, figsize=(6,9), gridspec_kw={'height_ratios': [4, 1]})
        plt.style.use("default")
        fig.suptitle('Center test', fontsize=12)
        plt.sca(ax1)
        s_g2, events, gates = gate_elipse(s_g1,
                                          channels=["SSC-A", "FSC-A"],
                                          center=center,
                                          a=0.3,
                                          b=0.2,
                                          t=30,
                                          events=events,
                                          gates=gates,
                                          plot=True,
                                          contour=True,
                                          plot_original=True,
                                          log=True)
        plt.sca(tb1)
        table(events, gates)
        plt.tight_layout()
        plt.savefig(output_path+'centertest.png', dpi=300)
        plt.close()
        print('Check file in input folder')
        done = input('Satisfied? y/n: ')
        if done != 'y':
            while True:
                try:
                    cx = float(input('new x coordinate: '))
                    break
                except KeyboardInterrupt:
                    sys.exit(1)
                except:
                    print('Center coordinate must be numerical value')
            
            while True:
                try:
                    cy = float(input('new y coordinate: '))
                    break
                except KeyboardInterrupt:
                    sys.exit(1)
                except:
                    print('Center coordinate must be numerical value')
            center = (cx, cy)
    os.remove(output_path+'centertest.png')
    return(center, channels)

def myExpFunc(x, a, b):
    return a * np.power(x, b)

def func(df, plot=False, plot_start=-1):
    
    xdata = np.array(df["FL1-A"].values.tolist())
    ydata = np.array(df["FL3-A"].values.tolist())
    popt, pcov = curve_fit(myExpFunc, xdata, ydata, maxfev=5000)
    
    residuals = ydata - myExpFunc(xdata, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r2 = 1-(ss_res/ss_tot)
    if plot:
        text = (f"y = {round(popt[0],2)}*x^{round(popt[1],2)}" + '\n' +
                "r^2 = " + str(round(r2, 2)))
        newX = np.logspace(plot_start, 2.9, base=10)
        plt.plot(newX, myExpFunc(newX, *popt), color='k', linestyle='--', marker='')  
        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                         lw=0, alpha=0)] * 2
        plt.legend(handles, [text], loc='best', fontsize='medium', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0)
    
    return(popt, r2)

def read_enhance_fcs(filename):
    s_raw = FlowCal.io.FCSData(filename)
    s = FlowCal.transform.to_rfi(s_raw)
    s_e = FlowCal.gate.high_low(s, channels=['SSC-A', 'FSC-A'])
    events = [len(s_e)]
    gates = ["No gate"]
    return(s, s_e, events, gates)

def gate_density(s, channels, gate_fraction, events, gates, plot=False, plot_mode='scatter', xscale='logicle', yscale='logicle'):
    gate_number=len(gates)
    gate_name="Density gate ("+str(gate_number)+")"
    gates.append(gate_name)
    s_g = FlowCal.gate.density2d(s,
                                 channels=channels,
                                 gate_fraction=gate_fraction)
    s_events = len(s_g)
    events.append(s_events)
    if plot:
        FlowCal.plot.density2d(s_g,
                          channels=channels,
                          mode=plot_mode,
                          xscale=xscale,
                          yscale=yscale)
    return(s_g, events, gates)

def gate_elipse(s, channels, center, a, b, t, events, gates, plot=False, plot_original=False, plot_mode='scatter', xscale='logicle', yscale='logicle', contour=False, log=False):
    gate_number=len(gates)
    gate_name="Eliptical gate ("+str(gate_number)+")"
    gates.append(gate_name)    
    s_g, m_g, contour_g = FlowCal.gate.ellipse(s,
                                               channels=channels,
                                               log=log,
                                               center=center,
                                               a=a,
                                               b=b,
                                               theta=t/180.*np.pi,
                                               full_output=True)
    s_events = len(s_g)
    events.append(s_events)
    if plot_original:
        s_plot = s
    else:
        s_plot = s_g    
    if plot:
        FlowCal.plot.density2d(s_plot,
                          channels=channels,
                          mode=plot_mode,
                          xscale=xscale,
                          yscale=yscale)
    if contour:
        for g in contour_g:
            plt.plot(g[:,0], g[:,1], color='k', linewidth=1.25)         
        
    return(s_g, events, gates)

def gate_low(s, gate_channel, low, events, gates, plot=False, plot_original=False, channels=0, plot_mode='scatter', xscale='logicle', yscale='logicle', contour=False):
    gate_number=len(gates)
    gate_name="Linear gate ("+str(gate_number)+")"
    gates.append(gate_name)     
    s_g = FlowCal.gate.high_low(s, 
                                channels=gate_channel, 
                                low=low)
    s_events = len(s_g)
    events.append(s_events)
    if plot_original:
        s_plot = s
    else:
        s_plot = s_g
        
    if plot:
        FlowCal.plot.density2d(s_plot,
                          channels=channels,
                          mode=plot_mode,
                          xscale=xscale,
                          yscale=yscale)
    if contour:
        plt.axvline(low, color="k", linewidth=1.25)
    return(s_g, events, gates) 

def gate_high(s, gate_channel, high, events, gates, plot=False, plot_original=False, channels=0, plot_mode='scatter', xscale='logicle', yscale='logicle', contour=False):
    gate_number=len(gates)
    gate_name="Linear gate ("+str(gate_number)+")"
    gates.append(gate_name)     
    s_g = FlowCal.gate.high_low(s, 
                                channels=gate_channel, 
                                high=high)
    s_events = len(s_g)
    events.append(s_events)
    if plot_original:
        s_plot = s
    else:
        s_plot = s_g
   
    if plot:
        FlowCal.plot.density2d(s_plot,
                          channels=channels,
                          mode=plot_mode,
                          xscale=xscale,
                          yscale=yscale)
    if contour:
        plt.axvline(high, color="k", linewidth=1.25)
    return(s_g, events, gates)

def adaptive_ellipse(s, channel1, channel2):
    df_s = pd.DataFrame(data=s, columns=s.channels)
    low_gate = df_s[channel1].mean()*1.2
    
    popt, r2 = popt, r2 = func(df_s)
    
    s_g = FlowCal.gate.high_low(s, channels=channel1, low=low_gate)
    df = pd.DataFrame(data=s_g, columns=s.channels)
    values, bins = np.histogram(df[channel1], bins=np.logspace(start=np.log10(1), stop=np.log10(10000000), num=50))
    order = np.argsort(values)[::-1]
    center_fl1 = bins[order[0]]
    center_fl3 = popt[0]*center_fl1+popt[1]
    center_fl3 = popt[0]*center_fl1**popt[1]
    values, bins = np.histogram(df[channel2], bins=np.logspace(start=np.log10(1), stop=np.log10(10000000), num=50))
    order = np.argsort(values)[::-1]
    center_fl3 = bins[order[0]]
    center_fl3 = df_s[channel2].mean()
    coordinate_fl1 = center_fl1**10
    coordinate_fl3 = popt[0]*coordinate_fl1+popt[1]
    coordinate_fl3 = popt[0]*coordinate_fl1**popt[1]
    x_fl1 = coordinate_fl1-center_fl1
    y_fl3 = coordinate_fl3-center_fl3
    log_x_center = math.log(center_fl1, 10)
    log_y_center = math.log(center_fl3, 10)
    
    log_x = math.log(x_fl1, 10)
    log_y = math.log(y_fl3, 10)
    angle = math.atan2(log_y, log_x)
    center = (log_x_center, log_y_center)
    
    return(center, angle)

def plot_only(s, channels, plot_mode='scatter', xscale='logicle', yscale='logicle'):
    FlowCal.plot.density2d(s,
                          channels=channels,
                          mode=plot_mode,
                          xscale=xscale,
                          yscale=yscale) 

def table(events, gates, fontsize=14):
    plt.axis("off")
    columns = ("Events", "Percent gated", "From previous")
    rows = gates
    data = [[events[0], "100.0%", "100.0%"]]
    for i in range(1,len(events)):
        data.append([events[i], str(round(100*events[i]/events[0], 2)) + "%", str(round(100*events[i]/events[i-1], 2)) + "%"])
        
    plt.table(cellText=data,
              rowLabels=rows,
              colLabels=columns,
              loc="center",
              cellLoc='center',
              colWidths = [0.25, 0.25, 0.25,0.25],
              fontsize=fontsize   
              )

def count_means(df, results, channel1, channel2, name):
    ch1mean = df[channel1].mean()
    ch2mean = df[channel2].mean()
    
    ch2mean_div_ch1mean = ch2mean/ch1mean
    
    divs = [ch2/ch1 for ch2, ch1 in zip(df[channel2], df[channel1])]
    
    divs_mean = statistics.mean(divs)
    divs_stdev = statistics.stdev(divs)
    
    results[name + " mean " + channel2].append(ch2mean)
    results[name + " mean " + channel1].append(ch1mean)
    results[name + " mean " + channel2 + " div mean " + channel1].append(ch2mean_div_ch1mean)
    results[name + " mean " + channel2 + " div mean " + channel1 + " one-by-one"].append(divs_mean)
    results[name + " mean " + channel2 + " div mean " + channel1 + " one-by-one stdev"].append(divs_stdev)

    return(results) 

def plot_files(file, results, center, channel1, channel2, high_gate, adapt_ellipse, fc_plot_path, sup_path):
    high_gate=float(high_gate)
    fig, ((ax1, ax2), (tb1, tb2), (ax3, ax4), (tb3, tb4)) = plt.subplots(4,2, figsize=(12,18), gridspec_kw={'height_ratios': [4, 1, 4, 1]})
    plt.style.use("default")
    filename = file.split("/")[-1].split(".")[0]
    fig.suptitle(filename, y=0.91, fontsize=25)
    results["file"].append(filename)
    
    s, s_g1, events, gates = read_enhance_fcs(file)
    plt.sca(ax1)    
    s_g2, events, gates = gate_elipse(s_g1,
                                       channels=["SSC-A", "FSC-A"],
                                       center=center,
                                       a=0.4,
                                       b=0.2,
                                       t=30,
                                       events=events,
                                       gates=gates,
                                       plot=True,
                                       contour=True,
                                       plot_original=True,
                                       log=True)     
    plt.sca(tb1)
    table(events, gates)
    
    plt.sca(ax2)
    plot_only(s_g2,
              channels = ["SSC-A", "FSC-A"],
              )

           
    plt.sca(tb2)
    
    s_g4 = s_g2
    table(events, gates)
            
    plt.sca(ax3)
    s_g5, events, gates = gate_low(s_g4,
                                   gate_channel=channel1,
                                   low=high_gate,
                                   events=events,
                                   gates=gates,
                                   plot=True,
                                   channels=[channel1, channel2],
                                   contour=True,
                                   plot_original=True
                                   )    

    plt.sca(tb3)
    table(events, gates)
    

    plt.sca(ax4)
    #FL1 FL3 ELLIPSE DENSITY TEST
    if adapt_ellipse==True:
        center_fl, angle = adaptive_ellipse(s_g5, channel1, channel2)
        s_g5, events, gates = gate_elipse(s_g5,
                                          channels=[channel1, channel2],
                                          center=center_fl,
                                          a=1,
                                          b=0.25,
                                          t=angle*180./np.pi,
                                          plot=True,
                                          plot_original=True,
                                          contour=True,
                                          log=True,
                                          events=events,
                                          gates=gates)
    else:
        plot_only(s_g5,
                  channels=[channel1, channel2])
    df_g5 = pd.DataFrame(data=s_g5, columns=s.channels)
    results = count_means(df_g5, results, channel1, channel2, "gated")
    results["gated events"].append(events[-1])
    

    plt.sca(tb4)
    table(events, gates)

    fc_savename = fc_plot_path + file.split(".")[0].split('/')[-1] + ".png"
    fig.savefig(fc_savename, dpi=300)  
    plt.close(fig)
    
    sup_savename = sup_path + file.split(".")[0].split('/')[-1] + ".png"
    sup_fig = plt.figure(figsize=(8, 4))
    sup_fig.suptitle(filename, fontsize=12)
    ax1 = plt.subplot2grid((13,20),(0, 0), rowspan=9, colspan=8)
    s_p1, m_p1, contour_p1 = FlowCal.gate.ellipse(s_g1,
                                                  ["SSC-A", "FSC-A"],
                                                  log=True,
                                                  center=center,
                                                  a=0.4,
                                                  b=0.2,
                                                  theta=30/180.*np.pi,
                                                  full_output=True)
    FlowCal.plot.density2d(s_g1,
                           channels=["SSC-A", "FSC-A"],
                           mode='scatter',
                           xscale='logicle',
                           yscale='logicle')
    for i in contour_p1:
        plt.plot(i[:,0], i[:,1], color='k', linewidth=1.25) 
    ax2 = plt.subplot2grid((13,20),(0, 11), rowspan=9, colspan=8)
    plot_only(s_g4,
              channels=['FL1-A', 'FL3-A'])
    plt.axvline(high_gate, color="k", linewidth=1.25)
    if adapt_ellipse:
        s_p2, m_p2, contour_p2 = FlowCal.gate.ellipse(s_g1,
                                                      [channel1, channel2],
                                                      log=True,
                                                      center=center_fl,
                                                      a=1,
                                                      b=0.25,
                                                      theta=angle,
                                                      full_output=True)
        for i in contour_p2:
            plt.plot(i[:,0], i[:,1], color='k', linewidth=1.25) 
    ax3 = plt.subplot2grid((13,20),(11,6), rowspan=9, colspan=8)
    table(events, gates, fontsize=20)    
    sup_fig.savefig(sup_savename, dpi=150)
    
    plt.close(sup_fig)
    return(results)

def analysis_iterator(input_dir, analyte, files, center, channel1, channel2, high_gate, adapt_ellipse):
    fc_plot_path, sup_path, results_path = make_output(analyte)
    results = results_format(channel1, channel2)
    pbar = tqdm(range(len(files)))
    for i in pbar:
        pbar.set_description("Processing %s" % files[i].split('/')[-1])
        sleep(0.1)
        file = files[i]
        results = plot_files(file, results, center, channel1, channel2, high_gate, adapt_ellipse, fc_plot_path, sup_path)
    savename = results_path + analyte + ".csv" 
    df = pd.DataFrame(results)
    df.to_csv(savename, index=False)
    
def main(args):
    print('')
    argument_check(args)
    input_dir = args.Input
    analyte = args.Analyte
    files = make_filelist(input_dir, analyte)
    print('Input directory:', input_dir)
    print('Analyte:', analyte)
    print('')
    
    center, channel1, channel2, high_gate, adapt_ellipse = parameters(args, input_dir, analyte, files)
    
    print('')
    print('____________________________')
    print('Running')
    print('SSC/FSC center:', center)
    print('Fluorescent channel for expression:', channel1)
    print('Fluorescent channel for binding:', channel2)
    print('Expression gate cutoff:', high_gate)
    if adapt_ellipse:
        print('Using the experimental feature to gate out the lower population')
    print('')
    
    analysis_iterator(input_dir, analyte, files, center, channel1, channel2, high_gate, adapt_ellipse)
    
    print('')
    print('Done!')
    print('')

if __name__ == '__main__':
    prog = "Epitope mapping initial analysis program"
    description = """Version 0.2, author Maximilian Karlander \n
    This script is designed for the initial analysis of FCS-files from the CytoFlex instrument for epitopemapping
    """
    epilog = "I did good, didn't I?"
    parser = argparse.ArgumentParser(prog=prog, 
                                     description=description, 
                                     epilog=epilog,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-d", "--Defaults", help="Define new defaults", action='store_true')
    parser.add_argument("-i", "--Input", help="Input directory", nargs='?', const="data/", type=str, default="data/")
    parser.add_argument("-a", "--Analyte", help="Analyte or antibody")
    parser.add_argument("-ch1", "--Channel1", help="Fluorescent channel for expression")
    parser.add_argument("-ch2", "--Channel2", help="Fluorescent channel for binding")
    parser.add_argument("-hi", "--HighGate", help="Gate to remove non expressing population")
    parser.add_argument("-dp", "--Double", help="!Experimental! Use this to gate out the lower population of two", action='store_true')
    args = parser.parse_args()  
    main(args)