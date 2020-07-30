"""
Created on Mon Jun 29 14:33:06 2020

@author: flavio martins

CSC method

"""

# Header ######################################################################

import numpy  as np
import pickle
import matplotlib.pyplot as plt

###############################################################################

def readLine(line):
    Value = []
    for word in line.split():
        try:
            Value.append(float(word))
        except ValueError:
            pass
        
    return Value
            
###############################################################################          

def readConfigFile(Flow):  
    # Name of the .txt file:
    name_of_file = './input/' + Flow + "_post_process_setup.txt"  
    
    fd = open(name_of_file,'r')  
    Lines = fd.readlines()
    for line in Lines:
        if 'dy' in line:
            dy = readLine(line)
        elif 'dz' in line:
            dz = readLine(line)
        elif 'Axis' in line:
            Axis = readLine(line)
        elif 'max_n_cases' in line:
            max_n_cases = readLine(line)
    
    return dy, dz, Axis, max_n_cases


###############################################################################

Flow = input('Name the flow (swirling_jet or ahmed_body): ')

# Read file setup:
mySetup = readConfigFile(Flow)

# Read setup file with configurations:
dy   = mySetup[0][0]
dz   = mySetup[1][0]
Axis = mySetup[2]
max_n_cases = int(mySetup[3][0])
    
# Read case from txt file:
file_with_cases = open('./output/' + Flow + '_cases' + '.txt', 'r')



count = 0
# For all lines in the text file containing cases:
for line in file_with_cases:
    # Text file containing test cases:
    file_name = './output/solution_' + Flow + line.strip() + '.pkl'
        
    # Process data for a given Np and Nt:
    with open(file_name, 'rb') as handle:
        
        # x,y - positions at t=0
        Data = pickle.load(handle)
        y = Data[1][0]
        z = Data[2][0]
        
        # Color and size of the scattered data:
        c = Data[6]*np.ones(np.shape(y))

        # Create figure:
        fig, ax = plt.subplots(figsize=(5.65/2.54, 5.65/2.54))
        sp = ax.scatter(y, z, c=c, s=1.5, cmap='jet')
        # fig.colorbar(sp) # add colorbar
        
        # use LaTeX fonts in the plot
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        
        # Set the axis:
        plt.xlabel('$y/D$', fontsize=10)
        plt.ylabel('$z/D$', fontsize=10)
        ax.set(xlim=(Axis[0], Axis[1]), ylim=(Axis[2], Axis[3]))
        plt.xticks(np.arange(Axis[0], 1.1*Axis[1], step=dy), fontsize=11)
        plt.yticks(np.arange(Axis[2], 1.1*Axis[3], step=dz), fontsize=11)
    
        # Plot results:     
        fig.savefig('./output/plot_' + Flow + '_' + line.strip() + '.pdf', bbox_inches='tight')
        
        plt.show()

        # Prevent the code from plotting all figures:
        count += 1
        if count >= max_n_cases:
            break
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    