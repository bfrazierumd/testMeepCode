#circuit radiation pattern testing
from __future__ import division
import math
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import h5py

from circuitFunctions import createCircuitGeometry

#constants 
c = 299792458

#work in cm (set scaling factor a = 1 cm for normalization)
a = 1e-2;

#setup the baseline parameters
fIn = 1e9;
fWidth = .5e9;
resolution = 30; #number of points per normalization length

#normalization
wavelength = c/fIn;
f = a/wavelength
nPointsPerWavelength = resolution/a*wavelength;
df = fWidth*a/c;
tScale = 2*resolution;

#get the geometry and physical layout for the simulation 
geometry,cell,pml_layers, p1Loc, nfXPos, nfYPos, nfXSize, nfYSize, board, buffer = createCircuitGeometry(a, resolution);

bx = board.size.x + board.center.x;
by = board.size.y + board.center.y;

print("Board X: " , bx, "Board Y: ", by)

maxPathLength = np.maximum(nfXSize.norm(),nfYSize.norm());
t = maxPathLength/c;
tau = t*c/a;
nTimeSteps = tScale*tau;

#output normalized values
print("Wavelength: ", wavelength, "m")
print("Normalized Frequency: ", f)
print("Normalized Frequency Width: ", df)
print("Number of Points Per Wavelength: ", nPointsPerWavelength)
print("Stop Time (seconds): ", t)
print("Number of Time Steps: ", nTimeSteps)
print("Normalized Stop Time: ", tau)


#setup the source
src_cmpt = mp.Ez
sources = [mp.Source(src=mp.GaussianSource(f,fwidth=df),
                     center=p1Loc,
                     component=src_cmpt)]

if src_cmpt == mp.Ex:
    symmetries = [mp.Mirror(mp.X,phase=-1),
                  mp.Mirror(mp.Y,phase=+1)]
    src_string="|Ex|";
    src_output = mp.output_efield_x;
elif src_cmpt == mp.Ey:
    symmetries = [mp.Mirror(mp.X,phase=+1),
                  mp.Mirror(mp.Y,phase=-1)]
    src_string="|Ey|";
    src_output = mp.output_efield_y;
elif src_cmpt == mp.Ez:
    symmetries = [mp.Mirror(mp.X,phase=+1),
                  mp.Mirror(mp.Y,phase=+1)]
    src_string="|Ez|";
    src_output = mp.output_efield_z;
	
#setup the overall simulation            
sim = mp.Simulation(cell_size=cell,
                    resolution=resolution,
                    sources=sources,
                    symmetries=symmetries,
                    boundary_layers=pml_layers,
                    geometry=geometry) 
                    
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez", mp.at_every(0.6, src_output)),
        until=tau)
        
eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Dielectric Constant")
plt.show()

x = np.linspace(0, bx*resolution, 100)
b = np.ones(100)*resolution;
o = np.ones(100);
y = np.linspace(0, by*resolution, 100)

e_data = sim.get_array(center=mp.Vector3(), size=cell, component=src_cmpt)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(e_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)

plt.plot(x + buffer*b,buffer*b + by*b ,color = 'black');
plt.plot(x + buffer*b,buffer*b, color = 'black');
plt.plot(buffer*b, y + buffer*b, color = 'black');
plt.plot(buffer*b + bx*b, y + buffer*b, color = 'black');
plt.xlabel("X")
plt.ylabel("Y")
plt.colorbar();
plt.title(src_string)
plt.show()