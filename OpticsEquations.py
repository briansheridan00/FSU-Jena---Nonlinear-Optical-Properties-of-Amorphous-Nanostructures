### Fresnel and Fabry-Perot equations ###


#Author: Brian Sheridan
#Credit: Navid Daryakar, Christin David
#Date: 25/June/2021


### What is this file? ###
#
#Two main functions are included in this file,
#the Fresnel equations for a single optical
#interface, and the Fabry-Perot equations for
#two optical interfaces. These functions
#take in material and optical parameters
#and output the reflection and transmission
#coefficients, and also the reflectance and
#transmittance (and absorptance).
#
#Material parameters may be taken from a related
#Python file 'MaterialParameters.py', which
#contains a database of the optical parameters
#of certain materials and can calculate the
#relative permittivity using the Drude-Lorentz
#model.
#
#See the MaterialParameters file for further details.



import numpy as np
import matplotlib.pyplot as plt


### Un-comment the following lines and
### insert the directory path of the
### MaterialParameters file, if it is
### not located in the same folder as
### this file.

#import sys
#sys.path.inset(0, '/Users/you/YourFolderHere')
import MaterialParameters as mp



e = 1.602176565e-19; # Electron charge [A s] or [J].
h = 6.62606957e-34;  # Planck's constant in [J s}.
c = 299792458e9;     # Speed of light in [nm / s}.
nprism = 1           # Refractive index of prism.



#Function to convert angular frequencies [eV] to wavelengths [nm].
def LfromW(frequency):
    lam = (h * c) / (frequency * e)
    return lam

#Function to convert wavelengths [nm] to angular frequencies [eV].
def WfromL(wavelength):
    W = (h * c) / (wavelength * e)
    return W



#Define Fresnel Equation function (used for one interface).
def Fresnel(polarisation="p",    #P-polarisation default.
            angle=0,             #Enter angle in degrees. (If array -> Must have same length as frequency/wavelength array.)
            wavelength=None,
            frequency=None,      #Input either wavelengths or frequencies.
            e0 = 1,              #Relative permittivity left hand side material (default air).
            e1 = -5.365+2.25j,   #Rel. Perm. of RHS material (Default Gold at w=3).
            extract="R"):        #Parameter the function returns (default Reflectance).

    a = angle * (np.pi) / 180    #Convert degrees to radians.

    if wavelength is None and frequency is not None:       #If the frequencies have been inputted.
        wavelength = LfromW(frequency)

    elif frequency is None and wavelength is not None:      #If the wavelengths have been inputted.
        wavelength = wavelength

    else:                        #If neither wavelengths nor frequencies are inputted.
        return print("ERROR: Please enter valid wavelengths or frequencies.")

    kx  = (2 * np.pi / wavelength) * nprism * np.sin(a)      #X projection of k-vector.
    kz0 = np.sqrt( e0 * (2*np.pi/wavelength)**2 - (kx)**2 )  #LHS medium - perpendicular projection.
    kz1 = np.sqrt( e1 * (2*np.pi/wavelength)**2 - (kx)**2 )  #RHS medium - perpendicular projection.


    #In order to conserve energy, the imaginary part
    #of the perpendicular wavenumber must not be negative.
    #We check if the imaginary part is positive in each
    #case and apply a sign change if negative. This is done
    #for each of the kz0 and kz1 wavenumbers, while
    #also checking whether it is a complex number, an
    #array of complex numbers or an array of arrays of complex
    #numbers. This is done in the cases of returning a single
    #reflectance/transmittance value, an array of values
    #dependent on a single variable (1D graph) or a
    #meshgrid array used to plot the dependence on 2 variables
    #(meshgrid arrays). These enhance the plotting flexibility.

    if isinstance(kz0, int) or isinstance(kz0, complex):
        if np.imag(kz0) < 0:
            kz0 = -kz0

    elif isinstance(kz0, np.ndarray) or isinstance(kz0, list):
        for i in range(0,len(kz0)):
            if isinstance(kz0[i], int) or isinstance(kz0[i], complex):
                if np.imag(kz0[i]) < 0:
                    kz0[i] = -kz0[i]
            elif isinstance(kz0[i], np.ndarray) or isinstance(kz0[i], list):
                for j in range(0, len(kz0[i])):
                    if np.imag(kz0[i][j]) < 0:
                        kz0[i][j] = -kz0[i][j]


    if isinstance(kz1, int) or isinstance(kz1, complex):
        if np.imag(kz1) < 0:
            kz1 = -kz1

    elif isinstance(kz1, np.ndarray) or isinstance(kz1, list):
        for i in range(0,len(kz1)):
            if isinstance(kz1[i], int) or isinstance(kz1[i], complex):
                if np.imag(kz1[i]) < 0:
                    kz1[i] = -kz1[i]
            elif isinstance(kz1[i], np.ndarray) or isinstance(kz1[i], list):
                for j in range(0, len(kz1[i])):
                    if np.imag(kz1[i][j]) < 0:
                        kz1[i][j] = -kz1[i][j]


    if polarisation=="p":
        r01 = ( (kz0/e0) - (kz1/e1) ) / ( (kz0/e0) + (kz1/e1) )  #r coefficient.
        t01 = ( 2 * np.sqrt(e0*e1) * kz0 ) / ( e1*kz0 + e0*kz1 ) #t coefficient.

    elif polarisation=="s":
        r01 = ( (kz0) - (kz1) ) / ( (kz0) + (kz1) )  #r coefficient.
        t01 = ( 2 * kz0 ) / ( kz0 + kz1 )            #t coefficient.

    R = r01 * np.conjugate(r01)  #Reflectance.
    R = np.real(R)               #Discard the imaginary part (which is 0 anyway).
    T = (1/np.cos(a)) * np.sqrt((e1/e0) - np.sin(a)**2) * t01 * np.conjugate(t01) #Transmittance.
    T = np.real(T)               #Discard the imaginary part (which is 0 anyway).

    if extract == "r":    #Return r coefficient if chosen.
        return r01
    elif extract == "t":  #Return t coefficient if chosen.
        return t01
    elif extract == "R":  #Return Reflectance if chosen.
        return R
    elif extract == "T":  #Return Transmittance if chosen.
        return T



#Define Fabry-Perot Equation function (used for two interfaces).
def FabryPerot(polarisation="p", #P-polarisation default.
            angle=0,             #Enter angle in degrees. Must have same length as frequency/wavelength array.
            wavelength=None,
            frequency=None,      #Input either wavelengths or frequencies.
            e0 = 1,              #Relative permittivity left hand side material (default air).
            e1 = -5.365+2.25j,   #Rel. Perm. of middle material (Default Gold at w=3).
            e2 = (1.54)**2,      #Rel.Perm of RHS material (Default SiO2).
            thickness = 200,     #Thickness of slab in nanometres.
            extract="R"):        #Parameter the function returns (default Reflectance).

    a = angle * (np.pi) / 180    #Convert degrees to radians.

    if wavelength is None and frequency is not None:   #If the frequencies have been inputted.
        wavelength = LfromW(frequency)

    elif frequency is None and wavelength is not None: #If the wavelengths have been inputted.
        wavelength = wavelength

    else:                        #If neither wavelengths nor frequencies are inputted.
        return print("ERROR: Please enter valid wavelengths or frequencies.")

    kx  = (2 * np.pi / wavelength) * nprism * np.sin(a)      #X projection of k-vector.
    kz0 = np.sqrt( e0 * (2*np.pi/wavelength)**2 - (kx)**2 )  #LHS medium - perpendicular projection.
    kz1 = np.sqrt( e1 * (2*np.pi/wavelength)**2 - (kx)**2 )  #Middle medium - perpendicular projection.
    kz2 = np.sqrt( e2 * (2*np.pi/wavelength)**2 - (kx)**2 )  #RHS medium - perpendicular projection.



    #In order to conserve energy, the imaginary part
    #of the perpendicular wavenumber must not be negative.
    #We check if the imaginary part is positive in each
    #case and apply a sign change if negative. This is done
    #for each of the kz0, kz1 and kz2 wavenumbers, while
    #also checking whether it is a complex number, an
    #array of complex numbers or an array of arrays of complex
    #numbers. This is done in the cases of returning a single
    #reflectance/transmittance value, an array of values
    #dependent on a single variable (1D graph) or a
    #meshgrid array used to plot the dependence on 2 variables
    #(meshgrid arrays). These enhance the plotting flexibility.

    if isinstance(kz0, int) or isinstance(kz0, complex):
        if np.imag(kz0) < 0:
            kz0 = -kz0

    elif isinstance(kz0, np.ndarray) or isinstance(kz0, list):
        for i in range(0,len(kz0)):
            if isinstance(kz0[i], int) or isinstance(kz0[i], complex):
                if np.imag(kz0[i]) < 0:
                    kz0[i] = -kz0[i]
            elif isinstance(kz0[i], np.ndarray) or isinstance(kz0[i], list):
                for j in range(0, len(kz0[i])):
                    if np.imag(kz0[i][j]) < 0:
                        kz0[i][j] = -kz0[i][j]


    if isinstance(kz1, int) or isinstance(kz1, complex):
        if np.imag(kz1) < 0:
            kz1 = -kz1

    elif isinstance(kz1, np.ndarray) or isinstance(kz1, list):
        for i in range(0,len(kz1)):
            if isinstance(kz1[i], int) or isinstance(kz1[i], complex):
                if np.imag(kz1[i]) < 0:
                    kz1[i] = -kz1[i]
            elif isinstance(kz1[i], np.ndarray) or isinstance(kz1[i], list):
                for j in range(0, len(kz1[i])):
                    if np.imag(kz1[i][j]) < 0:
                        kz1[i][j] = -kz1[i][j]


    if isinstance(kz2, int) or isinstance(kz2, complex):
        if np.imag(kz2) < 0:
            kz2 = -kz2

    elif isinstance(kz2, np.ndarray) or isinstance(kz2, list):
        for i in range(0,len(kz2)):
            if isinstance(kz2[i], int) or isinstance(kz2[i], complex):
                if np.imag(kz2[i]) < 0:
                    kz2[i] = -kz2[i]
            elif isinstance(kz2[i], np.ndarray) or isinstance(kz2[i], list):
                for j in range(0, len(kz2[i])):
                    if np.imag(kz2[i][j]) < 0:
                        kz2[i][j] = -kz2[i][j]



    if polarisation=="p":
        r01 = ( (kz0/e0) - (kz1/e1) ) / ( (kz0/e0) + (kz1/e1) )  #r coefficient - first interface.
        r12 = ( (kz1/e1) - (kz2/e2) ) / ( (kz1/e1) + (kz2/e2) )  #r coefficient - second interface.
        t01 = ( 2 * np.sqrt(e0*e1) * kz0 ) / ( e1*kz0 + e0*kz1 ) #t coefficient - first interface.
        t12 = ( 2 * np.sqrt(e1*e2) * kz1 ) / ( e2*kz1 + e1*kz2 ) #t coefficient - second interface.

    elif polarisation=="s":
        r01 = ( (kz0) - (kz1) ) / ( (kz0) + (kz1) )  #r coefficient - First interface.
        r12 = ( (kz1) - (kz2) ) / ( (kz1) + (kz2) )  #r coefficient - Second interface.
        t01 = ( 2 * kz0 ) / ( kz0 + kz1 )            #t coefficient - First interface.
        t12 = ( 2 * kz1 ) / ( kz1 + kz2 )            #t coefficient - Second interface.

    sval = 2 * kz1 * thickness
    pval = np.exp(sval * 1j)
    qval = np.exp(sval * 1j / 2)

    r012 = (r01 + r12 * pval ) / (1 + r01 * r12 * pval)
    t012 = (t01 * t12 * qval ) / (1 + r01 * r12 * pval)

    R = r012 * np.conjugate(r012)  #Reflectance.
    R = np.real(R)                 #Discard the imaginary part (which is 0 anyway).

    T = (1/np.cos(a)) * np.sqrt((e2/e0) - np.sin(a)**2) * t012 * np.conjugate(t012) #Transmittance.
    T = np.real(T)                 #Discard the imaginary part (which is 0 anyway).

    if extract == "r":    #Return r coefficient if chosen.
        return r012
    elif extract == "t":  #Return t coefficient if chosen.
        return t012
    elif extract == "R":  #Return Reflectance if chosen.
        return R
    elif extract == "T":  #Return Transmittance if chosen.
        return T
    elif extract == "A":  #Return Absorptance if chosen.
        return 1 - R - T
