### Material Parameters Database ###


#Author: Brian Sheridan
#Credit: Christin David, Navid Daryakar
#Date: 25/June/2021


### What is this file? ###
#
#This file is a database in Python of material parameters used
#to investigate optical transmission through
#different materials. This file stores relevant material
#parameters and can calculate the relative permittivity of
#a material using the Drude-Lorentz model, and can then plot
#the real and imaginary parts (using the numpy and
#matplotlib modules).


### How does this file work? ###
#
#A Python class for a general material is made,
#in which the relevent parameters are set to zero by default.
#Each material is an instance of this Python class
#for which the relevant material parameters
#(which have been found experimentally) are updated.
#Once the material is instantiated, you may call on the
#material parameters/functions in the form '[material].[parameter]'.
#Included in the Material class is a function to calculate the
#relative permittivity, and a function to plot the relative permittivity.
#The relative permittivity of each material is calculated using
#the Drude-Lorentz model and experimentally found material parameters,
#from S.G. Rodrigo et. al.
#(https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.075401)


### How to run this file. ###
#
#This materials database is written in the programming language
#Python, and therefore needs Python to work.
#Python may be downloaded from the official website
#(https://www.python.org/) which comes with a
#development environment called IDLE. You can open this
#document in IDLE to run it.
#
#Warning: if using the plotting function, this file
#uses the modules numpy and matplotlib. To install these
#(on Windows), open the command prompt once Python has been
#installed, then enter 'pip3 install numpy' and
#likewise for matplotlib.
#
#To run this file directly from the command prompt,
#change the working directory to that where the file is saved
#by entering 'cd Users/you/YourFolder', then run the file by
#entering 'python3 MaterialParameters.py'
#
#Alternatively you may use the Jupyter notebook environment
#which comes with relevant modules and packages pre-installed.
#
#In order to use this file in another Python script,
#make use of the code below. Note that if this file is
#in the same folder as your code, only the last line is needed.
#
#import sys
#sys.path.insert(0, 'Users/you/YourFolder/')
#import MaterialParameters
#
#Once Python is downloaded, the sys package is imported in which
#the directory of this Materials Database file may be linked.
#Then the Materials Database can be imported and used in other files.




#Create Materials class.
class Material:
    #Set the material parameters to zero by default.
    #The relevant parameters are updated for each material when the class is instantiated.
    def __init__(self, er=0, deleps1=0, deleps2=0, \
                 deleps3=0, deleps4=0, wp=0, wp1=0,\
                 wp2=0, wp3=0, gam0=0, gam1=0, gam2=0,\
                 gam3=0, Omega1=0, Omega2=0, Omega3=0,\
                 Omega4=0, Gamma1=0, Gamma2=0, Gamma3=0,\
                 Gamma4=0):

        self.er = er            #Relative permittivity.
        self.deleps1 = deleps1  #Delta of relative permittivity.
        self.deleps2 = deleps2
        self.deleps3 = deleps3
        self.deleps4 = deleps4
        self.wp = wp            #Plasma frequency.
        self.wp1 = wp1
        self.wp2 = wp2
        self.wp3 = wp3
        self.gam0 = gam0        #Damping constant.
        self.gam1 = gam1
        self.gam2 = gam2
        self.gam3 = gam3
        self.Omega1 = Omega1    #Oscillator frequency.
        self.Omega2 = Omega2
        self.Omega3 = Omega3
        self.Omega4 = Omega4
        self.Gamma1 = Gamma1    #Permittivity Relation.
        self.Gamma2 = Gamma2
        self.Gamma3 = Gamma3
        self.Gamma4 = Gamma4

    #Permittivity calculated using the Drude Lorentz model (Rodrigo et. al.)
    def epsilon(self, freq):
        eps = (self.er) \
        - (self.wp1)**2 / (freq *(freq + self.gam1*1j)) \
        - (self.wp2)**2 / (freq *(freq + self.gam2*1j)) \
        - (self.wp3)**2 / (freq *(freq + self.gam3*1j)) \
        - (self.deleps1)*(self.Omega1**2) /(freq**2 - self.Omega1**2 + freq*self.Gamma1*1j)\
        - (self.deleps2)*(self.Omega2**2) /(freq**2 - self.Omega2**2 + freq*self.Gamma2*1j)\
        - (self.deleps3)*(self.Omega3**2) /(freq**2 - self.Omega3**2 + freq*self.Gamma3*1j)\
        - (self.deleps4)*(self.Omega4**2) /(freq**2 - self.Omega4**2 + freq*self.Gamma4*1j)
        return eps

    #Plotting function for the real and imaginary parts of the permittivity.
    def plot(self, start=1.5, end=40, Npoints=1000, title="Relative Permittivity"):
        import numpy as np
        w = np.linspace(start, end, Npoints)
        import matplotlib.pyplot as plt
        Perm = self.epsilon(freq = w)
        plt.plot(w, np.real(Perm), label="Real")
        plt.plot(w, np.imag(Perm), label="Imaginary")
        plt.legend(loc = "best")
        plt.title(title)
        plt.xlabel("Angular Frequency (w)")
        plt.ylabel("Relative Permittivity")
        plt.grid(True)


#Create Material - Gold.
Au = Material()
Au.er=5.967;      #Relative Permittivity.
Au.wp1=8.729;     #First Plasma Frequency. Units [eV]
Au.gam0=0.065;    #First Damping constant. Units [eV]
Au.gam1=0.065;    #Second Damping Constant. Units [eV]
Au.deleps1=1.09;  #Delta epsilon.
Au.Omega1=2.684;  #Oscillator frequency. Units [eV]
Au.Gamma1=0.433;  # Units [eV]
#Au.diff = 1.9**2 * Au.wp1; # after Ashcroft & Mermin in nm^2 eV
#Au.diff = 1.3**2 * Au.wp1; # after Blader et al. in nm^2 eV


#Create Material - Silver.
Ag = Material()
Ag.er=4.6;
Ag.wp1=9.0;
Ag.gam0=0.07;
Ag.gam1=0.07;
Ag.deleps1=1.10;
Ag.Omega1=4.9;
Ag.Gamma1=1.20;
#Ag.diff = 1.9**2 * Ag.wp1; # after Ashcroft & Mermin in nm^2 eV
#Ag.diff = 1.1**2 * Ag.wp1; # after Blader et al. in nm^2 eV


#Create Material - Copper.
Cu = Material()
Cu.er=1.0;
Cu.wp1=8.212;
Cu.gam1=0.030;
Cu.deleps1=84.49;
Cu.Omega1=0.291;
Cu.Gamma1=0.378;
Cu.deleps2=1.395;
Cu.Omega2=2.957;
Cu.Gamma2=1.056;
Cu.deleps3=3.018;
Cu.Omega3=5.300;
Cu.Gamma3=3.213;
Cu.deleps4=0.598;
Cu.Omega4=11.18;
Cu.Gamma4=4.305;


#Create Material - Aluminium.
Al = Material()
Al.er=1.0;
Al.wp1=10.83;
Al.gam0=0.047;
Al.gam1=0.047;
Al.deleps1=1940.97;
Al.Omega1=0.162;
Al.Gamma1=0.333;
Al.deleps2=4.706;
Al.Omega2=1.544;
Al.Gamma2=0.312;
Al.deleps3=11.390;
Al.Omega3=1.808;
Al.Gamma3=1.351;
Al.deleps4=0.558;
Al.Omega4=3.473;
Al.Gamma4=3.382;


#Create Material - Nickel.
Ni = Material()
Ni.er=1.0;
Ni.wp1=4.621;
Ni.gam0=0.021;
Ni.gam1=0.021;
Ni.wp2=6.929;
Ni.gam2=1.771;
Ni.wp3=7.062;
Ni.gam3=3.443;
Ni.deleps1=2.1;
Ni.Omega1=1.458;
Ni.Gamma1=1.021;
Ni.deleps2=1.2;
Ni.Omega2=3.443;
Ni.Gamma2=2.410;


#Create Material - Chromium.
Cr = Material()
Cr.er=1.0;
Cr.wp1=4.406;
Cr.gam0=0.047;
Cr.gam1=0.047;
Cr.deleps1=1191.85;
Cr.Omega1=0.121;
Cr.Gamma1=3.175;
Cr.deleps2=  58.79;
Cr.Omega2=0.543;
Cr.Gamma2=1.305;
Cr.deleps3=  34.21;
Cr.Omega3=1.970;
Cr.Gamma3=2.676;
Cr.deleps4=  1.238;
Cr.Omega4=8.775;
Cr.Gamma4=1.335;


#Create Material - Tungsten.
W = Material()
W.er=1.0;
W.wp1=5.955;
W.gam0=0.027;
W.gam1=0.027;
W.wp2=2.286;
W.gam2=0.335;
W.deleps1=12.0;
W.Omega1=0.984;
W.Gamma1=0.590;
W.deleps2=14.4;
W.Omega2=2.066;
W.Gamma2=1.653;
W.deleps3=12.9;
W.Omega3=4.132;
W.Gamma3=2.479;


#Create Material - Vanadium Nitride.
VN = Material()
VN.er=1.47;
VN.wp1=8.04;
VN.gam0=0.91;
VN.gam1=0.91;
VN.deleps1=8.75;
VN.Omega1=6.63;
VN.Gamma1=9.50;


#Create Material - Molybdenum Nitride.
MoN = Material()
MoN.er=1.75;
MoN.wp1=9.69;
MoN.gam0=2.76;
MoN.gam1=2.76;
MoN.deleps1=5.49;
MoN.Omega1=5.18;
MoN.Gamma1=2.1;


#Create Material - Tungsten Nitride (1/2 Data Sets).
WN1 = Material()
WN1.er=1.00;
WN1.wp1=10.49;
WN1.gam0=3.19;
WN1.gam1=3.19;
WN1.deleps1=7.78;
WN1.Omega1=5.40;
WN1.Gamma1=2.45;


#Create Material - Tungsten Nitride (2/2 Data Sets).
WN2 = Material()
WN2.er=1.75;
WN2.wp1=8.15;
WN2.gam0=6.80;
WN2.gam1=6.80;
WN2.deleps1=16.95;
WN2.Omega1=3.78;
WN2.Gamma1=13.37;
WN2.deleps2= 0.83;
WN2.Omega2=5.85;
WN2.Gamma2=1.92;


#Create Material - Titanium Nitride (1/14 Data Sets)
TiN1 = Material()
TiN1.er=2.35;
TiN1.wp1=7.25;
TiN1.gam0=0.64;
TiN1.gam1=0.64;
TiN1.deleps1=0.23;
TiN1.Omega1=2.07;
TiN1.Gamma1=0.54;
TiN1.deleps2=4.04;
TiN1.Omega2=5.65;
TiN1.Gamma2=3.62;


#Create Material - Titanium Nitride (2/14 Data Sets)
TiN2 = Material()
TiN2.er=8.78;
TiN2.wp1=9.77;
TiN2.gam0=0.35;
TiN2.gam1=0.35;
TiN2.deleps1=0.40;
TiN2.Omega1=2.20;
TiN2.Gamma1=0.82;
TiN2.deleps2=5.28;
TiN2.Omega2=8.23;
TiN2.Gamma2=4.04;


#Create Material - Titanium Nitride (3/14 Data Sets)
TiN3 = Material()
TiN3.er=3.00;
TiN3.wp1=8.08;
TiN3.gam0=0.86;
TiN3.gam1=0.86;
TiN3.deleps1=3.54;
TiN3.Omega1=5.14;
TiN3.Gamma1=3.04;


#Create Material - Titanium Nitride (4/14 Data Sets)
TiN4 = Material()
TiN4.er=3.84;
TiN4.wp1=6.42;
TiN4.gam0=0.86;
TiN4.gam1=0.86;
TiN4.deleps1=3.12;
TiN4.Omega1=3.77;
TiN4.Gamma1=1.55;


#Create Material - Titanium Nitride (5/14 Data Sets)
TiN5 = Material()
TiN5.er=2.13;
TiN5.wp1=7.05;
TiN5.gam0=0.63;
TiN5.gam1=0.63;
TiN5.deleps1=0.38;
TiN5.Omega1=3.71;
TiN5.Gamma1=1.33;
TiN5.deleps2=4.35;
TiN5.Omega2=5.88;
TiN5.Gamma2=4.28;


#Create Material - Titanium Nitride (6/14 Data Sets)
TiN6 = Material()
TiN6.er=3.87;
TiN6.wp1=7.21;
TiN6.gam0=0.81;
TiN6.gam1=0.81;
TiN6.deleps1=2.69;
TiN6.Omega1=4.48;
TiN6.Gamma1=2.64;
TiN6.deleps2=1.27;
TiN6.Omega2=5.36;
TiN6.Gamma2=0.88;


#Create Material - Titanium Nitride (7/14 Data Sets)
TiN7 = Material()
TiN7.er=2.67;
TiN7.wp1=5.71;
TiN7.gam0=0.17;
TiN7.gam1=0.17;
TiN7.deleps1=0.20;
TiN7.Omega1=2.26;
TiN7.Gamma1=0.72;
TiN7.deleps2=2.31;
TiN7.Omega2=4.89;
TiN7.Gamma2=3.88;


#Create Material - Titanium Nitride (8/14 Data Sets)
TiN8 = Material()
TiN8.er=1.85;
TiN8.wp1=7.75;
TiN8.gam0=0.35;
TiN8.gam1=0.35;
TiN8.deleps1=0.10;
TiN8.Omega1=3.68;
TiN8.Gamma1=0.71;
TiN8.deleps2=5.62;
TiN8.Omega2=6.26;
TiN8.Gamma2=3.88;


#Create Material - Titanium Nitride (9/14 Data Sets)
TiN9 = Material()
TiN9.er=1.87;
TiN9.wp1=6.93;
TiN9.gam0=0.59;
TiN9.gam1=0.59;
TiN9.deleps1=0.18;
TiN9.Omega1=3.69;
TiN9.Gamma1=0.94;
TiN9.deleps2=4.88;
TiN9.Omega2=5.97;
TiN9.Gamma2=4.88;


#Create Material - Titanium Nitride (10/14 Data Sets)
TiN10 = Material()
TiN10.er=1.20;
TiN10.wp1=4.49;
TiN10.gam0=1.38;
TiN10.gam1=1.38;
TiN10.deleps1=0.55;
TiN10.Omega1=3.76;
TiN10.Gamma1=2.08;
TiN10.deleps2=1.64;
TiN10.Omega2=6.67;
TiN10.Gamma2=5.55;


#Create Material - Titanium Nitride (11/14 Data Sets)
TiN11 = Material()
TiN11.er=1.20;
TiN11.wp1=6.94;
TiN11.gam0=0.73;
TiN11.gam1=0.73;
TiN11.deleps1=3.96;
TiN11.Omega1=5.52;
TiN11.Gamma1=4.12;
TiN11.deleps2=1.04;
TiN11.Omega2=7.82;
TiN11.Gamma2=2.20;


#Create Material - Titanium Nitride (12/14 Data Sets)
TiN12 = Material()
TiN12.er=1.96;
TiN12.wp1=6.36;
TiN12.gam0=0.74;
TiN12.gam1=0.74;
TiN12.deleps1=0.15;
TiN12.Omega1=3.48;
TiN12.Gamma1=0.76;
TiN12.deleps2=7.01;
TiN12.Omega2=5.79;
TiN12.Gamma2=5.93;


#Create Material - Titanium Nitride (13/14 Data Sets)
TiN13 = Material()
TiN13.er=3.18;
TiN13.wp1=8.05;
TiN13.gam0=0.95;
TiN13.gam1=0.95;
TiN13.deleps1=1.47;
TiN13.Omega1=3.99;
TiN13.Gamma1=2.54;
TiN13.deleps2=2.51;
TiN13.Omega2=5.34;
TiN13.Gamma2=2.26;


#Create Material - Titanium Nitride (14/14 Data Sets)
TiN14 = Material()
TiN14.er=3.24;
TiN14.wp1=7.07;
TiN14.gam0=0.44;
TiN14.gam1=0.44;
TiN14.deleps1=4.22;
TiN14.Omega1=5.08;
TiN14.Gamma1=3.35;


#Create Material - Zirconium Nitride (1/7 Data Sets).
ZrN1 = Material()
ZrN1.er=1.38;
ZrN1.wp1=8.08;
ZrN1.gam0=0.41;
ZrN1.gam1=0.41;
ZrN1.deleps1=5.11;
ZrN1.Omega1=7.98;
ZrN1.Gamma1=5.91;


#Create Material - Zirconium Nitride (2/7 Data Sets).
ZrN2 = Material()
ZrN2.er=1.79;
ZrN2.wp1=7.11;
ZrN2.gam0=0.51;
ZrN2.gam1=0.51;
ZrN2.deleps1=3.41;
ZrN2.Omega1=6.29;
ZrN2.Gamma1=3.15;


#Create Material - Zirconium Nitride (3/7 Data Sets).
ZrN3 = Material()
ZrN3.er=2.99;
ZrN3.wp1=8.07;
ZrN3.gam0=0.18;
ZrN3.gam1=0.18;
ZrN3.deleps1=2.57;
ZrN3.Omega1=5.83;
ZrN3.Gamma1=1.54;


#Create Material - Zirconium Nitride (4/7 Data Sets).
ZrN4 = Material()
ZrN4.er=3.22;
ZrN4.wp1=7.35;
ZrN4.gam0=0.40;
ZrN4.gam1=0.40;
ZrN4.deleps1=2.34;
ZrN4.Omega1=6.51;
ZrN4.Gamma1=2.73;


#Create Material - Zirconium Nitride (5/7 Data Sets).
ZrN5 = Material()
ZrN5.er=5.41;
ZrN5.wp1=8.16;
ZrN5.gam0=0.47;
ZrN5.gam1=0.47;
ZrN5.deleps1=0.82;
ZrN5.Omega1=4.52;
ZrN5.Gamma1=1.25;


#Create Material - Zirconium Nitride (6/7 Data Sets).
ZrN6 = Material()
ZrN6.er=2.44;
ZrN6.wp1=7.02;
ZrN6.gam0=0.73;
ZrN6.gam1=0.73;
ZrN6.deleps1=3.97;
ZrN6.Omega1=6.33;
ZrN6.Gamma1=4.34;


#Create Material - Zirconium Nitride (7/7 Data Sets).
ZrN7 = Material()
ZrN7.er=2.90;
ZrN7.wp1=7.66;
ZrN7.gam0=0.33;
ZrN7.gam1=0.33;
ZrN7.deleps1=1.92;
ZrN7.Omega1=5.27;
ZrN7.Gamma1=1.51;


#Create Material - Hafnium Nitride (1/5 Data Sets).
HfN1 = Material()
HfN1.er=1.74;
HfN1.wp1=8.09;
HfN1.gam0=0.92;
HfN1.gam1=0.92;
HfN1.deleps1=3.41;
HfN1.Omega1=6.48;
HfN1.Gamma1=3.32;


#Create Material - Hafnium Nitride (2/5 Data Sets).
HfN2 = Material()
HfN2.er=1.30;
HfN2.wp1=7.30;
HfN2.gam0=0.12;
HfN2.gam1=0.12;
HfN2.deleps1=3.97;
HfN2.Omega1=4.44;
HfN2.Gamma1=2.62;


#Create Material - Hafnium Nitride (3/5 Data Sets).
HfN3 = Material()
HfN3.er=2.09;
HfN3.wp1=5.33;
HfN3.gam0=0.53;
HfN3.gam1=0.53;
HfN3.deleps1=0.35;
HfN3.Omega1=3.96;
HfN3.Gamma1=1.40;


#Create Material - Hafnium Nitride (4/5 Data Sets).
HfN4 = Material()
HfN4.er=4.32;
HfN4.wp1=5.20;
HfN4.gam0=0.47;
HfN4.gam1=0.47;
HfN4.deleps1=61.25;
HfN4.Omega1=0.70;
HfN4.Gamma1=0.86;
HfN4.deleps2=1.83;
HfN4.Omega2=4.77;
HfN4.Gamma2=1.52;


#Create Material - Hafnium Nitride (5/5 Data Sets).
HfN5 = Material()
HfN5.er=6.58;
HfN5.wp1=8.01;
HfN5.gam0=0.36;
HfN5.gam1=0.36;
HfN5.deleps1=0.04;
HfN5.Omega1=6.05;
HfN5.Gamma1=16.77;


#Create Material - Niobium Nitride (1/3 Data Sets).
NbN1 = Material()
NbN1.er=1.92;
NbN1.wp1=7.17;
NbN1.gam0=2.42;
NbN1.gam1=2.42;
NbN1.deleps1=1.47;
NbN1.Omega1=5.32;
NbN1.Gamma1=3.56;


#Create Material - Niobium Nitride (2/3 Data Sets).
NbN2 = Material()
NbN2.er=2.93;
NbN2.wp1=9.15;
NbN2.gam0=1.81;
NbN2.gam1=1.81;
NbN2.deleps1=2.22;
NbN2.Omega1=5.33;
NbN2.Gamma1=2.39;


#Create Material - Niobium Nitride (3/3 Data Sets).
NbN3 = Material()
NbN3.er=2.69;
NbN3.wp1=8.89;
NbN3.gam0=1.97;
NbN3.gam1=1.97;
NbN3.deleps1=2.74;
NbN3.Omega1=6.19;
NbN3.Gamma1=4.43;


#Create Material Tantalum Nitride (1/7 Data Sets).
TaN1 = Material()
TaN1.er=1.54;
TaN1.wp1=3.41;
TaN1.gam0=0.99;
TaN1.gam1=0.99;
TaN1.deleps1=26.70;
TaN1.Omega1=1.80;
TaN1.Gamma1=4.30;
TaN1.deleps2= 4.85;
TaN1.Omega2=7.50;
TaN1.Gamma2=7.60;


#Create Material Tantalum Nitride (2/7 Data Sets).
TaN2 = Material()
TaN2.er=3.49;
TaN2.wp1=2.88;
TaN2.gam0=1.12;
TaN2.gam1=1.12;
TaN2.deleps1=13.49;
TaN2.Omega1=1.63;
TaN2.Gamma1=1.65;
TaN2.deleps2= 3.88;
TaN2.Omega2=5.42;
TaN2.Gamma2=3.68;


#Create Material Tantalum Nitride (3/7 Data Sets).
TaN3 = Material()
TaN3.er=1.80;
TaN3.wp1=3.92;
TaN3.gam0=1.88;
TaN3.gam1=1.88;
TaN3.deleps1=8.96;
TaN3.Omega1=1.92;
TaN3.Gamma1=2.55;
TaN3.deleps2=4.63;
TaN3.Omega2=7.42;
TaN3.Gamma2=6.38;


#Create Material Tantalum Nitride (4/7 Data Sets).
TaN4 = Material()
TaN4.er=2.32;
TaN4.wp1=5.55;
TaN4.gam0=0.80;
TaN4.gam1=0.80;
TaN4.deleps1=31.22;
TaN4.Omega1=1.45;
TaN4.Gamma1=2.25;
TaN4.deleps2= 3.73;
TaN4.Omega2=6.47;
TaN4.Gamma2=4.54;


#Create Material Tantalum Nitride (5/7 Data Sets).
TaN5 = Material()
TaN5.er=1.53;
TaN5.wp1=3.29;
TaN5.gam0=0.31;
TaN5.gam1=0.31;
TaN5.deleps1=42.69;
TaN5.Omega1=1.55;
TaN5.Gamma1=3.23;
TaN5.deleps2= 4.23;
TaN5.Omega2=6.27;
TaN5.Gamma2=2.28;


#Create Material Tantalum Nitride (6/7 Data Sets).
TaN6 = Material()
TaN6.er=2.10;
TaN6.wp1=9.45;
TaN6.gam0=2.85;
TaN6.gam1=2.85;
TaN6.deleps1=3.65;
TaN6.Omega1=7.24;
TaN6.Gamma1=5.61;


#Create Material Tantalum Nitride (7/7 Data Sets).
TaN7 = Material()
TaN7.er=1.68;
TaN7.wp1=3.48;
TaN7.gam0=1.02;
TaN7.gam1=1.02;
TaN7.deleps1=25.85;
TaN7.Omega1=1.80;
TaN7.Gamma1=4.21;
TaN7.deleps2=4.78;
TaN7.Omega2=7.01;
TaN7.Gamma2=6.95;
