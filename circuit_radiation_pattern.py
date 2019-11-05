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
fIn = 10e9;
fWidth = 5e9;
resolution = 30; #number of points per normalization length

#normalization
wavelength = c/fIn;
f = a/wavelength
nPointsPerWavelength = resolution/a*wavelength;
df = fWidth*a/c;
tScale = 2*resolution;


#output normalized values
print("Wavelength: ", wavelength)
print("Normalized Frequency: ", f)
print("Normalized Frequency Width: ", df)
print("Number of Points Per Wavelength: ", nPointsPerWavelength)

#get the geometry and physical layout for the simulation 
geometry,cell,pml_layers, p1Loc, nfXPos, nfYPos, nfXSize, nfYSize = createCircuitGeometry(a, resolution);

#setup the source
src_cmpt = mp.Ex
sources = [mp.Source(src=mp.GaussianSource(f,fwidth=df),
                     center=p1Loc,
                     component=src_cmpt)]
                     
if src_cmpt == mp.Ex:
    symmetries = [mp.Mirror(mp.X,phase=-1),
                  mp.Mirror(mp.Y,phase=+1)]
elif src_cmpt == mp.Ey:
    symmetries = [mp.Mirror(mp.X,phase=+1),
                  mp.Mirror(mp.Y,phase=-1)]
elif src_cmpt == mp.Ez:
    symmetries = [mp.Mirror(mp.X,phase=+1),
                  mp.Mirror(mp.Y,phase=+1)]

#setup the overall simulation            
sim = mp.Simulation(cell_size=cell,
                    resolution=resolution,
                    sources=sources,
                    symmetries=symmetries,
                    boundary_layers=pml_layers,
                    geometry=geometry)  
                  
 
#setup the near field box 
nearfield_box = sim.add_near2far(f, 0, 1, 
mp.Near2FarRegion(mp.Vector3(y = nfYPos), size=nfXSize),
mp.Near2FarRegion(mp.Vector3(y = -nfYPos), size=nfXSize, weight=-1),
mp.Near2FarRegion(mp.Vector3(x = nfXPos),size = nfYSize),
mp.Near2FarRegion(mp.Vector3(x = -nfXPos), size=nfYSize, weight=-1))
                            
#add the flux box
flux_box = sim.add_flux(f, 0, 1,
mp.FluxRegion(mp.Vector3(y = nfYPos), size=nfXSize),
mp.FluxRegion(mp.Vector3(y = -nfYPos), size=nfXSize, weight=-1),
mp.FluxRegion(mp.Vector3(x = nfXPos), size=nfYSize),
mp.FluxRegion(mp.Vector3(x = -nfXPos), size=nfYSize, weight=-1))

#run the simulation
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt, mp.Vector3(), 1e-8))


print("Getting Flux ...")
near_flux = mp.get_fluxes(flux_box)[0]

r = 1000/f      # half side length of far-field square box OR radius of far-field circle
res_ff = 1       # resolution of far fields (points/μm)
far_flux_box = (nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=r), size=mp.Vector3(2*r)), res_ff)[0]
               - nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=-r), size=mp.Vector3(2*r)), res_ff)[0]
               + nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(r), size=mp.Vector3(y=2*r)), res_ff)[0]
              - nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(-r), size=mp.Vector3(y=2*r)), res_ff)[0])

npts = 1000         # number of points in [0,2*pi) range of angles
angles = 2*math.pi/npts*np.arange(npts)

print("Computing Farfield ...")
E = np.zeros((npts,3),dtype=np.complex128)
H = np.zeros((npts,3),dtype=np.complex128)
for n in range(npts):
    ff = sim.get_farfield(nearfield_box,
                          mp.Vector3(r*math.cos(angles[n]),
                                     r*math.sin(angles[n])))
    E[n,:] = [np.conj(ff[j]) for j in range(3)]
    H[n,:] = [ff[j+3] for j in range(3)]

Px = np.real(np.multiply(E[:,1],H[:,2])-np.multiply(E[:,2],H[:,1]))
Py = np.real(np.multiply(E[:,2],H[:,0])-np.multiply(E[:,0],H[:,2]))
Pr = np.sqrt(np.square(Px)+np.square(Py))

far_flux_circle = np.sum(Pr)*2*np.pi*r/len(Pr)

print("flux:, {:.6f}, {:.6f}, {:.6f}".format(near_flux,far_flux_box,far_flux_circle))

#save the response to an h5 file
if src_cmpt == mp.Ex:
    hf=h5py.File('antenna_response_ex.h5','w');
elif src_cmpt == mp.Ey:
    hf=h5py.File('antenna_response_ey.h5','w');
elif src_cmpt == mp.Ez:
    hf=h5py.File('antenna_response_ez.h5','w');

hf.create_dataset('antennaPower', data=Pr);
hf.create_dataset('angles', data=angles);
hf.close

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Dielectric Constant")
plt.show()

ax = plt.subplot(111, projection='polar')
ax.plot(angles,Pr/max(Pr),'b-')
ax.set_rmax(1)
ax.set_rticks([0,0.5,1])
ax.grid(True)
ax.set_rlabel_position(22)
plt.show()

