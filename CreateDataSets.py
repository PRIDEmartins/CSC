"""
Created on Mon Jun 29 14:33:06 2020

@author: flavio martins

CSC method

"""

# Header ######################################################################

import numpy  as np
import pickle
import random
import progressbar

from scipy.io import loadmat

###############################################################################

# Input: <line> (string)
# Output: <Value> (float)
# This function reads all numbers separate by space in a string and returns an
# array containing all numbers. Ex: 'a 1 2 3' -> [1, 2, 3]

def readLine(line):
    Value = []
    for word in line.split():
        try:
            Value.append(float(word))
        except ValueError:
            pass
        
    return Value
            
###############################################################################

# Input: <Flow> (string)
# Output: <values associated with keywords> (tulep)
# This function reads txt files, search for keywords, and extracts all numbers
# next to them. Ex: search for (a) in 'a, 10, 20' -> [10, 20]

def readConfigFile(Flow):  
    # Name of the .txt file:
    name_of_file = './input/' + Flow + "_particles_selection.txt"  
    
    fd = open(name_of_file,'r')  
    Lines = fd.readlines()
    for line in Lines:
        if 'Nsub' in line:
            Nsub = readLine(line)
        elif 'Domain' in line:
            Domain = readLine(line)
        elif 'D' in line:
            D = readLine(line)
        elif 'Shift' in line:
            Shift = readLine(line)
        elif 'Bo' in line:
            Bo = readLine(line)
        elif 'Nt' in line:
            Nt = readLine(line)
        elif 'Np' in line:
            Np = readLine(line)
    
    return Nsub, Domain, D, Shift, Bo, Nt, Np

###############################################################################

#   Create Xflowmaps(Nt,Np) containing the values of location (x,y,z) and 
# velocities (u,v,z) of Np particles of equal track-length of Nt.
#   Particles are selected randomly, and data is normalized utilizing the 
# characteristic length D of the flow.

def createFlowMaps():
    # Load parameters:
    Nsub = int(mySetup[0][0])
    D = mySetup[2][0] #characteristic length of the flow
    Shift = mySetup[3]
    Bo = mySetup[4]   # boundaries for the values of y: [ymin, ymax]
    
    
    # This is very particular to the way this dataset is written. Other flow
    # data will require a different way of extracting the x,y,z,...,ID data
    if Flow == 'swirling_jet':
        # Load data from .mat file:
        loadData = loadmat('./input/'+ Flow + '.mat')
        Data = loadData['Data']
        
        # Transform data into arrays:
        x =  np.array(Data[:,1])/D 
        y =  np.array(Data[:,0])/D
        z =  np.array(Data[:,2])/D
        u =  np.array(Data[:,4])
        v =  np.array(Data[:,3])
        w =  np.array(Data[:,5])
        ID = np.array(Data[:,6])
        t =  np.array(Data[:,7])
        
    elif Flow == 'ahmed_body':
        # Load data from mat file:
        loadData = loadmat('./input/'+ Flow + '.mat')
        Data = loadData['Data_STB_Light']
        
        # Position
        x = (Data['coord_abs'][0][0:Nsub] + Shift[0])/D
        y = (Data['coord_abs'][1][0:Nsub] + Shift[1])/D
        z = (Data['coord_abs'][2][0:Nsub] + Shift[2])/D
        
        # Velocity
        u = (Data['vel_abs'][0][0:Nsub] + Shift[0])/D
        v = (Data['vel_abs'][1][0:Nsub] + Shift[1])/D
        w = (Data['vel_abs'][2][0:Nsub] + Shift[2])/D
        
        t  =  Data['timestep'][0][0:Nsub]
        ID =  Data['trackID'][0][0:Nsub]
        
    # Determine start and end-points of the tracks:
    a, starte = np.unique(ID, return_index=True)
    ende = np.array(starte[1:len(starte)]-1)
    ende = np.append(ende, ID[len(ID)-1])
    
    # Determine the length of the tracks:
    lengthe = ende-starte+1
    
    # Determine the number of particles (Np) and time-steps (Nt):
    Np = NPmax
    Nt = Lmin
    
    # Create flow maps:
    Xflowmap = np.empty((Nt,1))
    Yflowmap = np.empty((Nt,1))
    Zflowmap = np.empty((Nt,1))
    Uflowmap = np.empty((Nt,1))
    Vflowmap = np.empty((Nt,1))
    Wflowmap = np.empty((Nt,1))
    
    # Start extracting the tracks:
    pp = 0 # index of particles that respect the criteria
    for p in range(0,len(lengthe)-1):
        
        # If the track has a length equal or greater than Lmin:
        if lengthe[p] >= Lmin:
            
            # If the particle was in a y-range of [ymin,ymax]:
            if y[starte[p]] >= Bo[0] and y[starte[p]] <= Bo[1]:
                inDomain = 1
            else:
                inDomain = 0

            if inDomain:
                pos1 = starte[p]
                pos2 = starte[p]+Lmin
                
                # Read tracks from these locations:
                loc = np.arange(start=pos1, stop=pos2, step=1)
                
                # Transform the array x[loc] into a vertial array and append
                # the Xflowmap with the array of size 'loc'
                Xflowmap = np.append(Xflowmap, x[loc].reshape(-1,1), axis=1)
                Yflowmap = np.append(Yflowmap, y[loc].reshape(-1,1), axis=1)
                Zflowmap = np.append(Zflowmap, z[loc].reshape(-1,1), axis=1)
                Uflowmap = np.append(Uflowmap, u[loc].reshape(-1,1), axis=1)
                Vflowmap = np.append(Vflowmap, v[loc].reshape(-1,1), axis=1)
                Wflowmap = np.append(Wflowmap, w[loc].reshape(-1,1), axis=1)
                
                # Number of particles collected:
                pp = pp+1
                
    # Update the number of time-steps and particles:
    Nt = Lmin
    Np = min(pp-1,NPmax)

    # If at least 10 particles were found during the process, save data:
    if Np >= 10:
        # create a list of numbers in between 1 and pp
        list_ = np.arange(start=1, stop=pp, step=1) 
        # Randomly select Np particles in the list list_:
        particles = random.sample(list(list_),Np) 

        # Fiz the dimensions of the data:
        Xflowmap = Xflowmap[0:Nt,particles]
        Yflowmap = Yflowmap[0:Nt,particles]
        Zflowmap = Zflowmap[0:Nt,particles]
        Uflowmap = Uflowmap[0:Nt,particles]
        Vflowmap = Vflowmap[0:Nt,particles]
        Wflowmap = Wflowmap[0:Nt,particles]
    
    return Xflowmap, Yflowmap, Zflowmap, Uflowmap, Vflowmap, Wflowmap, Np, Nt

###############################################################################

# Input: <Flow> (string) and <Data> (tuple)
# Saves data from the tuple Data using the name as a reference. Data is of
# the form:
# Xflowmap, Yflowmap, Zflowmap, Uflowmap, Vflowmap, Wflowmap, Np, Nt
  
def saveData(Flow, Data):  
    Np = Data[6]
    Nt = Data[7]
    if Np >= 10: # data[6] = number of particles
        name_of_file = './output/' + Flow + '_flowmaps_Nt' + str(Nt) + '_Np' + str(Np) + '.pkl'
        with open(name_of_file, 'wb') as f:
            pickle.dump([Data[0],Data[1],Data[2],Data[3],Data[4],Data[5],Data[6],Data[7]], f)

###############################################################################

# Create txt file with lines containing the string 'Np<Np>_Nt<Nt>'.
# This file informs the code CSC.py which cases to post-process

def caseToFile(Np,Nt):
    if Np >= 10:
        nameOfFile = './output/' + Flow + '_cases.txt'
        with open(nameOfFile, 'a+') as f:
            f.write('Nt' + str(Nt) + '_Np' + str(Np) + '\n')
            f.close

# Input #######################################################################


Flow = input('Name the flow (swirling_jet or ahmed_body): ')

# Read setup file:
mySetup = readConfigFile(Flow)

# Values of track length given by the user:
track_lengths = mySetup[5]
# Number of tracer particles determined by the user:
num_of_particles = mySetup[6]

# <Number of cases> =
#      <different number of track-lengths> x <different number of particles>
n_p = len(track_lengths)
n_t = len(num_of_particles)
n_cases = n_p*n_t

print('\nCreating flowmaps.\n')

# Progress bar setup and initialization:
bar = progressbar.ProgressBar(maxval=n_cases, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
bar.start()


count = 1
for NL in range(0,n_p):
    for NP in range(0,n_t):
        NPmax = int(num_of_particles[NP])
        Lmin  = int(track_lengths[NL])
        
        # Read mySetup and create flowmaps (X,Y,Z,U,V,W):
        Data = createFlowMaps()
        
        # Save flowmaps to .pkl data:
        saveData(Flow, Data)       
        
        # Create .txt file with case (Np,Nt):
        caseToFile(Data[6], Data[7])
        
        # Print progress:
        bar.update(count)
        count += 1

# Close progress bar:
bar.finish()
print('\nDone.\n')

###############################################################################