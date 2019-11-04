from __future__ import division

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

#work in cm (set scaling factor a = 1 cm for normalization)
a = 1e-2;

UseEigenModeSolver = True;

#setup the baseline parameters
fIn = 10e9;
fWidth = 5e9;
resolution = 10; #number of points per normalization length
nfreq = 1000;
pmlThickness = 0.1;
buffer = 10;

#setup the board dimensions
bwidth = 4*2.54;
bheight = 3*2.54;
bthick = mp.inf;#0.5;
bdielectric = 4.4;

#setup the trace
twidth = 2.5*2.54;
theight = 0.1*2.54;
tthick = mp.inf;#0.01;
tconductivity = 58.108e6*a; # M-Siemens/m converted to the normalization

#constants 
c = 299792458

#normalization
wavelength = c/fIn;
f = a/wavelength
nPointsPerWavelength = resolution/a*wavelength;
df = fWidth*a/c;
tScale = 2*resolution;

#initialize the computational cell
cellX = bwidth + 2*buffer;
cellY = bheight + 2*buffer;
cellZ = 0;#bthick + 2;
cell = mp.Vector3(cellX,cellY,cellZ);

maxPathLength = 2;
t = maxPathLength/c;
tau = t*c/a;
nTimeSteps = tScale*tau;

#output normalized values
print("Wavelength: ", wavelength)
print("Normalized Frequency: ", f)
print("Normalized Frequency Width: ", df)
print("Number of Points Per Wavelength: ", nPointsPerWavelength)
print("Stop Time (seconds): ", t)
print("Number of Time Steps: ", nTimeSteps)
print("Normalized Stop Time: ", tau)

planeSize = mp.Vector3(0,theight,bthick);
inputPortLoc = mp.Vector3(-twidth/2,1,0);
outputPortLoc = mp.Vector3(-twidth/2,-1,0);
#initialize the board
board = mp.Block(mp.Vector3(bwidth,bheight,bthick),
                     center=mp.Vector3(0,0,0),
                     material=mp.Medium(epsilon=bdielectric));
                     
trace1 = mp.Block(mp.Vector3(twidth,theight,tthick),
                   center=mp.Vector3(0,1,0),
                   material=mp.metal);  
                   
trace2 = mp.Block(mp.Vector3(twidth,theight,tthick),
                   center=mp.Vector3(0,-1,0),
                   material=mp.metal);                    


geometry = [board, trace1, trace2];

pml_layers = [mp.PML(pmlThickness)]


if UseEigenModeSolver:                   
	sources = [mp.EigenModeSource(src=mp.GaussianSource(frequency=f,fwidth=df),
                                  size=planeSize,
                                  center=inputPortLoc,
                                  eig_band=1,
                                  eig_parity=mp.EVEN_Y+mp.ODD_Z,
                                  eig_match_freq=True)];
else:
	sources = [mp.Source(mp.ContinuousSource(frequency=f),
                     component=mp.Ez,
                    center=inputPortLoc)]

                                  
sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

if UseEigenModeSolver:                    
	plane1 = mp.Block(planeSize,center=inputPortLoc)
	plane2 = mp.Block(planeSize,center=outputPortLoc)
	mode1 = sim.add_mode_monitor(f, df, nfreq, mp.ModeRegion(volume=plane1,direction=mp.X))
	mode2 = sim.add_mode_monitor(f, df, nfreq, mp.ModeRegion(volume=plane2,direction=mp.X))

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z)),
        until=tau)

if UseEigenModeSolver:        
	S1 = sim.get_eigenmode_coefficients(mode1, [1], eig_parity=mp.EVEN_Y+mp.ODD_Z)
	S2 = sim.get_eigenmode_coefficients(mode2, [1], eig_parity=mp.EVEN_Y+mp.ODD_Z)
	eig_freqs = mp.get_eigenmode_freqs(mode1)

	S11 = []
	S12 = []
	S21 = []
	S22 = []
	ff = []

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
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Dielectric Constant")
plt.show()

x = np.linspace(0, bwidth*resolution, 100)
b = np.ones(100)*resolution;
y = np.linspace(0, bheight*resolution, 100)

ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.plot(x + buffer*b,buffer*b + bheight*b ,color = 'black');
plt.plot(x + buffer*b,buffer*b, color = 'black');
plt.plot(buffer*b, y + buffer*b, color = 'black');
plt.plot(buffer*b + bwidth*b, y + buffer*b, color = 'black');
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Ez Component")
plt.show()
