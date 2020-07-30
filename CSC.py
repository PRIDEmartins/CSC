"""
Created on Mon Jun 29 14:33:06 2020

@author: flavio martins

CSC method

"""

# Header ######################################################################

import numpy  as np
import pickle
import progressbar
import scipy.linalg as la
import os

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
    name_of_file = './input/' + Flow + "_CSC_setup.txt"  
    
    fd = open(name_of_file,'r')  
    Lines = fd.readlines()
    for line in Lines:
        if 'max_n_cases' in line:
            max_n_cases = readLine(line)
        elif 'Fs' in line:
            Fs = readLine(line)
    
    return max_n_cases, Fs

###############################################################################

def computeAdjacencyMatrix(Data):

    # print('Computing adjacency matrix:\n')
    
    # Load data:
    Np = Data[6]
    Nt = Data[7]
    Fs = mySetup[1][0]
    Xflowmap = Data[0]
    Yflowmap = Data[1]
    Zflowmap = Data[2]
    
    # Check if temp dir exists. If not, create it:
    #if not os.path.exists('temp'):
    #    os.makedirs('temp')
    #else:
    #    os.rmdir("temp")
    #    os.makedirs('temp')
    
    # Compute the adjacency matrix r:
    r = np.zeros((Np,Np))
    rsum = np.zeros((Np,Np))
    
    # Start progress bar:
    bar = progressbar.ProgressBar(max_value = Nt)

    
    # Compute adjacency matrix:
    count = 1
    for t in range(0, Nt):
        for p1 in range(0, Np):
            xp1 = Xflowmap[t,p1]
            yp1 = Yflowmap[t,p1]
            zp1 = Zflowmap[t,p1]
            for p2 in range(0, Np):
                xp2 = Xflowmap[t,p2]
                yp2 = Yflowmap[t,p2]
                zp2 = Zflowmap[t,p2]
                r[p1,p2] = ( (xp1-xp2)**2 + (yp1-yp2)**2 + (zp1-zp2)**2 )**0.5
        R = r + r.transpose()
        rsum = rsum + R
        r = np.zeros((Np,Np))
        
        # Save matrix R for every time-step:
        with open('C:/temp/R' + str(t) + '.pkl', 'wb') as f:
            pickle.dump([R], f)
        
        # Update progress:
        bar.update(count)
        count += 1
    
    # Mean value:
    rmean = rsum/Nt
    
    
    # Compute the adjacency matrix:
    Tf = (1/Fs)*Nt
    s = np.zeros(rmean.shape)
    for t in range(0, Nt):
        # For every time-step, load the matrix R:
        with open('C:/temp/R' + str(t) + '.pkl' , 'rb') as f:
            R = pickle.load(f)
        s = s + np.square(rmean-R[0])
    
    # Compute A and replace main diagonal for zeroes:
    A = (1/Tf**0.5)*(np.sqrt(s)/rmean)     
    np.fill_diagonal(A,0)
    
    # Free-up even more space:
    #del rmean s t Tf R
    
    return A
        
###############################################################################
    
def computeDegreeMatrix(A):   
    Np = Data[6]
    #Nt = Data[7]
    D = np.zeros((Np,Np))
    
    for i in range(0,Np):
        D[i,i] = np.sum(A[i,:])
    
    # Graph Laplacian:
    L = D-A
    # Normalize L:
    Dinv = np.sqrt(np.linalg.matrix_power(D,-1)) # D^-0.5
    L = Dinv @ L @ Dinv   # normalize L
    
    # del D A
    
    return L
  
############################################################################### 
    
def computeCSCsolution(L):
    
    # By definition, the CSC is the eigvec associated with the maximum
    # eig-value, therefore:
    eigvals, eigvecs = la.eig(L)
    max_eigval_pos = np.argmax(eigvals)
    CSCflowmap = np.real(eigvecs[:,max_eigval_pos])
    
    #del eigvals eigvecs max_eigval_pos
    
    return CSCflowmap
    
############################################################################### 
    
def saveResults():
        file_name = './output/solution_' + Flow + line.strip() + '.pkl'
    
        with open(file_name, 'wb') as f:
            pickle.dump([Data[0],Data[1],Data[2],Data[3],Data[4],Data[5],CSC], f)

###############################################################################
    
Flow = input('Name the flow (swirling_jet or ahmed_body): ')

# Read file setup:
mySetup = readConfigFile(Flow)
max_n_cases = int(mySetup[0][0])
    
# Read case from txt file:
file_with_cases = open('./output/' + Flow + '_cases' + '.txt', 'r')

count = 1
for line in file_with_cases:
    # Text file containing test cases:
    file_name = './output/' + Flow + '_flowmaps_' + line.strip() + '.pkl'
        
    print('\n\nSolving ' + str(count) + ' of ' + str(max_n_cases) + ' cases:')
        
    # for every line of the text file, process the data:
    with open(file_name, 'rb') as handle:
        Data = pickle.load(handle)
        A = computeAdjacencyMatrix(Data)
        L = computeDegreeMatrix(A)
        CSC = computeCSCsolution(L)
            
        saveResults()
    
    if count >= max_n_cases:
        break

    count += 1
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


