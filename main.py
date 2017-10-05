import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

# important things to note: positive residues will start
# from ph - pk and gradually ascend to 0 (when the x axis,
# ph value reaches that of pk. It will stay that way for the
# rest of the ph values on the x axis. For negative residues,
# the energy will start at 0, and then stay that way until ph
# (the x axis) matches the value of pk. After this event, the
# energy will decrease constantly until it reaches its final
# value of pk - ph. (last ph, of course, is 14)

class pk_shifts_prereqs(object):
        def __init__(self):
                self.pksol_dict = dict(ARG = 12.5, HIS = 6.98, LYS = 10.4,
                                ASP = 4.75, GLU = 4.75, TYR = 10.2)
                self.res_charges = dict(POSITIVE = ['ARG', 'HIS', 'LYS'],
                                NEGATIVE = ['ASP', 'GLU', 'TYR'])
                self.x_axis = np.linspace(0,14,15)

        def gather_data(self, pk_filepath):
                if not(os.path.exists(pk_filepath)):
                        sys.exit('SPECIFIED PATH FOR pK.out NOT FOUND')
                raw_data = pd.read_csv(pk_filepath, sep = '\s+', header = None,
                                        usecols = [0,1], names = ['Res','pK'])
                raw_data.set_index('Res', inplace = True)
                # using pd to filter out residues for which we have default
                # pkSol values. These residues are in the self.res_charges
                # dictionary.         
                self.data = raw_data[ [raw_data.index.values[i][0:3]
                                in np.ravel(list(p.res_charges.values())) for i in
                                range(len(raw_data.index))] ]
                # filter out amino acids that have pK values that are >14
                # or that are <0.
                self.data = self.data[ [self.data['pK'].values[i] not in
                                ['<0.0','>14.0'] for i in
                                range(len(self.data['pK']))] ]
                self.positive_res = self.data[ [self.data.index.values[i][0:3]
                                in self.res_charges['POSITIVE'] for i in
                                range(len(self.data.index))] ]
                self.negative_res = self.data[ [self.data.index.values[i][0:3]
                                in self.res_charges['NEGATIVE'] for i in 
                                range(len(self.data.index))] ]

def run():
    '''Instantiates the class above and creates graphs based on the 
    information provided
    '''
    p = pk_shifts_prereqs()
    p.gather_data('pK.out')

    to_plot_master_positive = np.zeros(15) # Declaration of Pos. Agg. Values
    to_plot_master_negative = np.zeros(15) # Declaration of Neg. Agg. Values
    # Declaration of the subplot figure
    f, axarr = plt.subplots(2, 2, figsize = (15, 15))

    # Plot Positive Residues Individually
    # Declaration of variables to plot
    for i in range(p.positive_res.size):    
        to_plot = np.zeros(15)
        a_pos_ind = [p.pksol_dict[p.positive_res.index.values[i][0:3]] - p.x_axis[j] if p.x_axis[j] <= 
             p.pksol_dict[p.positive_res.index.values[i][0:3]] else 0 for j in range(len(p.x_axis))]
        b_pos_ind = [p.positive_res.values[i].astype(float)[0] - p.x_axis[j] if p.x_axis[j] <= 
             p.positive_res.values[i].astype(float)[0] else 0 for j in range(len(p.x_axis))]
        to_plot_pos_ind = np.array(a_pos_ind) - np.array(b_pos_ind)
        axarr[0, 0].plot(p.x_axis, to_plot_pos_ind) # Calling the plot function for the above type, FOR EACH ITERATION
        to_plot_master_positive = to_plot_master_positive + to_plot_pos_ind # Aggregating all values from each iteration       
    axarr[0, 0].set_title('Pos. Ind.')

    # Plot Negative Residues Individually
    # Declaration of variables to plot
    for i in range(p.negative_res.size):
        to_plot = np.zeros(15)
        a_neg_ind = [p.pksol_dict[p.negative_res.index.values[i][0:3]] - p.x_axis[j] if p.x_axis[j] >= 
             p.pksol_dict[p.negative_res.index.values[i][0:3]] else 0 for j in range(len(p.x_axis))]
        b_neg_ind = [p.negative_res.values[i].astype(float)[0] - p.x_axis[j] if p.x_axis[j] >= 
             p.negative_res.values[i].astype(float)[0] else 0 for j in range(len(p.x_axis))]
        to_plot_neg_ind = np.array(b_neg_ind) - np.array(a_neg_ind)
        axarr[0, 1].plot(p.x_axis, to_plot_neg_ind) # Calling the plot function for the above type, FOR EACH ITERATION
        to_plot_master_negative = to_plot_master_negative + to_plot_neg_ind # Aggregating all values from each iteration
    axarr[0, 1].set_title('Neg. Ind.')

    # Plot Positive Residues Aggregated        
    axarr[1, 0].plot(p.x_axis, to_plot_master_positive) # Calling the plot function for the aggregate
    axarr[1, 0].set_title('Pos. Agg.')

    # Plot Negative Residues Aggregated
    axarr[1, 1].plot(p.x_axis, to_plot_master_negative) # Calling the plot function for the aggregate
    axarr[1, 1].set_title('Neg. Agg.')       

    for i in [0, 1]:
        for j in [0, 1]:
            axarr[i, j].grid(color = 'black', linestyle = '-', linewidth = 0.15)
            axarr[i, j].set_xlabel('$pH$')
            axarr[i, j].set_ylabel('$E$')
    
    f.savefig('analysis.png')
    
    # Plotting AGGREGATE of ALL residues
    new_fig = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    to_plot_master_all = to_plot_master_positive + to_plot_master_negative
    plt.plot(p.x_axis, to_plot_master_all) 
    plt.grid(color = 'blacK', linestyle = '-', linewidth = 0.15)
    plt.title('ALL Residues Aggregate pK-Shift')
    plt.xlabel('$pH$')
    plt.ylabel('$E$')
    
    plt.show()
    new_fig.savefig('pk_shift.png')

if __name__ == '__main__':
    p = pk_shifts_prereqs()
    #p.positive_res)
    run()
    
    
