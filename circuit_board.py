from __future__ import division

import meep as mp
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from meep.materials import Pt


#work in cm (set scaling factor a = 1 cm for normalization)
a = 1e-2;

#setup the baseline parameters
fIn = 1e9;
fWidth = 0.5e9;
resolution = 10; #number of points per normalization length
nfreq = 500;
pmlThickness = 0.1;
buffer = 10;
maxPathLength = 10;

#setup the board dimensions
bwidth = 4*2.54;
bheight = 3*2.54;
bthick = mp.inf;#0.5;
bdielectric = 4.4;

#setup the trace
twidth = 2.5*2.54;
theight = 0.1*2.54;
tthick = mp.inf;#0.01;
tconductivity = 58.108e4; # Siemens/cm

#constants 
c = 299792458

#normalization
wavelength = c/fIn;
f = a/wavelength
nPointsPerWavelength = resolution/a*wavelength;
df = fWidth*a/c;
t = maxPathLength/c;
tau = t*c/a;
tScale = 2*resolution;
nTimeSteps = tScale*tau;

#output normalized values
print("Wavelength: ", wavelength)
print("Normalized Frequency: ", f)
print("Normalized Frequency Width: ", df)
print("Number of Points Per Wavelength: ", nPointsPerWavelength)
print("Stop Time (seconds): ", t)
print("Number of Time Steps: ", nTimeSteps)
print("Normalized Stop Time: ", tau)


#initialize the computational cell
cellX = bwidth + 2*buffer;
cellY = bheight + 2*buffer;
cellZ = 0;#bthick + 2;
cell = mp.Vector3(cellX,cellY,cellZ);

planeSize = mp.Vector3(0,theight,bthick);
#initialize the board
board = mp.Block(mp.Vector3(bwidth,bheight,bthick),
                     center=mp.Vector3(0,0,0),
                     material=mp.Medium(epsilon=bdielectric));
                     
trace = mp.Block(mp.Vector3(twidth,theight,tthick),
                   center=mp.Vector3(0,0,bthick/2),
                   material=mp.metal);  

sources = [mp.EigenModeSource(src=mp.GaussianSource(frequency=f,fwidth=df),
                                  size=planeSize,
                                  center=mp.Vector3(-twidth/2,0,0),
                                  eig_band=1,
                                  eig_parity=mp.EVEN_Y+mp.ODD_Z,
                                  eig_match_freq=True)];
                                  
geometry = [board, trace];

pml_layers = [mp.PML(pmlThickness)]
sim = mp.Simulation(resolution=resolution,
                    cell_size=cell,
                    boundary_layers=pml_layers,
                    sources=sources,
                    geometry=geometry);
 
 
plane1 = mp.Block(planeSize,center=mp.Vector3(-twidth/2,0,0))
plane2 = mp.Block(planeSize,center=mp.Vector3(twidth/2,0,0))
mode1 = sim.add_mode_monitor(f, df, nfreq, mp.ModeRegion(volume=plane1,direction=mp.X))
mode2 = sim.add_mode_monitor(f, df, nfreq, mp.ModeRegion(volume=plane2,direction=mp.X))

sim.run(until_after_sources=tau);

S11 = []
S12 = []
S21 = []
S22 = []
ff = []

S1 = sim.get_eigenmode_coefficients(mode1, [1], eig_parity=mp.EVEN_Y+mp.ODD_Z)
S2 = sim.get_eigenmode_coefficients(mode2, [1], eig_parity=mp.EVEN_Y+mp.ODD_Z)
eig_freqs = mp.get_eigenmode_freqs(mode1)

for i in range(0,nfreq):
	S11 = np.append(S11,S1.alpha[0,i,0]);
	S12 = np.append(S12,S2.alpha[0,i,0]);
	ff = np.append(ff,eig_freqs[i]*c/a*1e-9);

plt.figure()
plt.plot(ff,np.abs(S11));
plt.plot(ff,np.abs(S12));
plt.xlabel("Frequency (GHz)");
plt.ylabel("|S|");
plt.legend("$S_{11}$", "$S_{12}$");
plt.show();

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.show()  

ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()                