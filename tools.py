import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c

'''
@author: Lucivaldo Aguiar 
References:
https://optics.ansys.com/hc/en-us/articles/360042819293-Coupled-ring-resonator-filters
Integrated ring resonator, Dominik G. Rabus
Synthesis of a parallel-coupled ring-resonator filter, Andrea Melloni
'''
def FSR2graph (FSR):
    #this graphic is valid only and if only m1 - m2 = 1
    nm = 1e-9
    um = 1e-6
    FSR1 = np.linspace(0.1*nm,15*nm, 100)
    FSR2 = FSR/(np.abs(1-(FSR/FSR1)))
    plt.plot(FSR1/nm,FSR2/nm)
    plt.ylim(min(FSR1/nm), max(FSR2/nm))
    plt.xlim(min(FSR1/nm), max(FSR1/nm))
    plt.xlabel('FSR1 [nm]')
    plt.ylabel('FSR2 [nm]')
    plt.title(rf'$FSR_2$($FSR ={FSR/nm}, FSR_1$)')
    plt.grid()
    plt.show()
    
def ringData(x, N,ng, FSR):
    """
       Computes key parameters for a ring resonator based on input parameters.

       Inputs:
           x   : Relation FSR/B (Free Spectral Range to Bandwidth ratio)
           N   : Number of rings
           ng  : Group index of the waveguide
           FSR : Free Spectral Range

       Output:
           A dictionary containing calculated values for FSR, K, Q, g, B, and L.
       """
    dictionary ={
        'FSR' : [],
        'K' : [],
        'Q' : [],
        'g' : [],
        'B' : [],
        'L' : []
        }
    #determining FSRs and ring total length
    if N > 1:
        m1 = 3
        m2 = 2
        FSR1 = FSR/m1
        FSR2 = (m1/m2)*FSR1
        dictionary['FSR'].extend([FSR1,FSR2])
        #ring total length
        l1 = (1550e-9**2)/(ng*FSR1)
        l2 = (1550e-9**2)/(ng*FSR2)
        dictionary['L'].extend([l1,l2])
    else:
        dictionary['FSR'].append(FSR)
        l = (1550e-9**2)/(ng*FSR)
        dictionary['L'].append(l)
        
    for n in range(1,3):
        #Computing parameters based on Melloni's method
        g = np.sqrt(2) * np.sin((2*n-1)/(2*N) * np.pi)
        B = FSR / x
        Q = FSR/(B*g)
        K = ((np.pi**2)/(2*Q**2))*(np.sqrt(1 + (4*Q**2)/(np.pi**2)-1))

        dictionary['g'].append(g)
        dictionary['Q'].append(FSR/(B*g))    
        dictionary['K'].append(K)
        dictionary['B'].append(B)
    #Computing coupling coefficient for the middle ring
    dictionary['K'].append(np.sqrt(0.25)*(dictionary['K'][0])**2)
    return dictionary

