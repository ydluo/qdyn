"""
helper functions to explore parameter space 
"""
import numpy as np

# Calculate cohesive zone (process zone) distribution
def cohesive_zone(G,Dc,b,sigma):
    """
    G: shear modulus
    Dc: characteristic slip distance
    b: evolution parameter
    sigma: normal stress

    Lb needs to be several times larger than the mesh resolution to
    resolve the seismic rupture propagation (Day et al., 2005)
    """

    Lb = G*Dc/(b*sigma)
    return Lb

# Calculate nucleation zone distribution
def nucleation_zone(G,Dc,b,sigma,a,b)
    """
    G: shear modulus
    Dc: characteristic slip distance
    b: evolution parameter
    a: state parameter
    sigma: normal stress

    Seismogenic width (W) needs to be several times larger than
    L_inf to resolve seismic nucleation (Rubin and Ampuero, 2005)
    """
    
    L_inf = G*Dc*b/(np.pi*sigma*((b-a)**2))
    return L_inf

# Calculate Ruina-Rice number

def ru_number(a,b,sigma,W,G,Dc):
     """
    G: shear modulus
    Dc: characteristic slip distance
    b: evolution parameter
    a: state parameter
    sigma: normal stress


    Dietrich-Ruina-Rice number (Ru) is the ratio of the seismogenic zone size 
    to a characteristic rupture nucleation dimension 
    (other approach to explore relationship between L_inf and W)

    Ru<1: stable creep
    Ru>>1: from stable ruptures to chaotic sequences with increasing Ru
    (Catania, 2019; Barbot et al., 2019)
    """

    Ru = ((b-a)*sigma*W)/(G*Dc)
    return Ru

 # Calculate Rb number
 def rb_number(a,b)
    """
    a: state parameter
    b: evolution parameter

    >> Ru: complex sequences
    << Ru: slows-slip events
    """
    Rb = (b-a)/b

    return Rb