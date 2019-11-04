from __future__ import division

import meep as mp
import numpy as np
import numpy as np
import matplotlib.pyplot as plt


#setup the baseline parameters
fIn = 10e9;
fWidth = 4e9;
resolution = 10
nfreq = 500;
#for scaling, let a = 1 cm
a = 1e-2;

wallWidth = 1
portLength = wallWidth
portHeight = 3
buffer = 6
pmlThickness = 2

waveguideLength = 30;
waveguideHeight = 20;
maxPathLength = 10;

#constants
c = 299792458

#normalization
wavelength = c/fIn
f = a/wavelength
nPointsPerWavelength = resolution/a*wavelength;
df = fWidth*a/c
t = maxPathLength/c;
tau = t*c/a;
tScale = 2*resolution;
nTimeSteps = tScale*tau;

dielectricConstant = 1

#output normalized values
print("Normalized Frequency: ", f)
print("Normalized Frequency Width: ", df)
print("Number of Points Per Wavelength: ", nPointsPerWavelength)
print("Stop Time (seconds): ", t)
print("Number of Time Steps: ", nTimeSteps)
print("Normalized Stop Time: ", tau)

##  define the cavity geometry
#setup the computational cell

cellX = waveguideLength + 2*wallWidth + 2*buffer;
cellY = waveguideHeight + 2*wallWidth + 2*buffer;

cell = mp.Vector3(cellX,cellY,0)

portSize = mp.Vector3(portLength,portHeight,0);
p1Loc = mp.Vector3(-waveguideLength/2-wallWidth/2,0,0);
p2Loc = mp.Vector3(waveguideLength/2+wallWidth/2,0,0);


waveguide = mp.Block(mp.Vector3(waveguideLength,waveguideHeight,mp.inf),
                     center=mp.Vector3(0,0,0),
                     material=mp.Medium(epsilon=dielectricConstant))
                     
topwall = mp.Block(mp.Vector3(waveguideLength+2*wallWidth,wallWidth,mp.inf),
                   center=mp.Vector3(0,-waveguideHeight/2-wallWidth/2,0),
                   material=mp.metal);                   
bottomwall = mp.Block(mp.Vector3(waveguideLength+2*wallWidth,wallWidth,mp.inf),
                   center=mp.Vector3(0,waveguideHeight/2+wallWidth/2,0),
                   material=mp.metal);   

leftwall = mp.Block(mp.Vector3(wallWidth,waveguideHeight+2*wallWidth,mp.inf),
                   center=mp.Vector3(-waveguideLength/2-wallWidth/2,0,0),
                   material=mp.metal);  
                   
rightwall = mp.Block(mp.Vector3(wallWidth,waveguideHeight+2*wallWidth,mp.inf),
                   center=mp.Vector3(waveguideLength/2+wallWidth/2,0,0),
                   material=mp.metal);  
                   
obstacle = mp.Cylinder(radius=2,height=1,axis=mp.Vector3(0,0,1),center=mp.Vector3(0,0,0),
                       material=mp.metal);
                       
obstacle1 = mp.Ellipsoid(center=mp.Vector3(-5,4,0), size=mp.Vector3(8,2,2),
                   e1=mp.Vector3(1,1,0), e2=mp.Vector3(0,1,0),e3=mp.Vector3(0,0,1),
                   material=mp.metal)
                   
obstacle2 = mp.Ellipsoid(center=mp.Vector3(6,-6,0), size=mp.Vector3(8,2,2),
                   e1=mp.Vector3(1,.5,0), e2=mp.Vector3(0,1,0),e3=mp.Vector3(0,0,1),
                   material=mp.metal)
                                     
port1 = mp.Block(portSize,
                 center=p1Loc,
                 material=mp.Medium(epsilon=dielectricConstant))
                 
port2 = mp.Block(portSize,
                 center=p2Loc,
                 material=mp.Medium(epsilon=dielectricConstant))
                                                      

geometry = [waveguide, topwall, bottomwall, leftwall, rightwall, port1, port2, obstacle1, obstacle2]


planeSize = mp.Vector3(0,portHeight,0)

                                     
sources = [mp.EigenModeSource(src=mp.GaussianSource(frequency=f,fwidth=df),
                                  size=planeSize,
                                  center=p1Loc,
                                  eig_band=1,
                                  eig_parity=mp.EVEN_Y+mp.ODD_Z,
                                  eig_match_freq=True)]
                                                    

   
pml_layers = [mp.PML(pmlThickness)]
sim = mp.Simulation(resolution=resolution,
                    cell_size=cell,
                    boundary_layers=pml_layers,
                    sources=sources,
                    geometry=geometry)
                   

plane1 = mp.Block(planeSize,center=p1Loc)
plane2 = mp.Block(planeSize,center=p2Loc)
mode1 = sim.add_mode_monitor(f, df, nfreq, mp.ModeRegion(volume=plane1,direction=mp.X))
mode2 = sim.add_mode_monitor(f, df, nfreq, mp.ModeRegion(volume=plane2,direction=mp.X))

sim.run(until_after_sources=tau)

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

#sim.reset_meep()
 
#print("Now running for Port 2") 

#sources = [mp.EigenModeSource(src=mp.GaussianSource(frequency=f,fwidth=df),
#                                  size=portSize,
 #                                 center=p2Loc,
  #                                eig_band=1,
   #                               eig_parity=mp.EVEN_Y+mp.ODD_Z,
    #                              eig_match_freq=True)]
  
#sim = mp.Simulation(resolution=resolution,
 #                   cell_size=cell,
  #                  boundary_layers=pml_layers,
   #                 sources=sources,
    #                geometry=geometry)
                                                    
#sim.run(until_after_sources=tau)


#S1 = sim.get_eigenmode_coefficients(mode1, [1], eig_parity=mp.EVEN_Y+mp.ODD_Z)
#S2 = sim.get_eigenmode_coefficients(mode2, [1], eig_parity=mp.EVEN_Y+mp.ODD_Z)

#for i in range(0,nfreq):
	#S21 = np.append(S21,S1.alpha[0,i,1]);
	#S22 = np.append(S22,S2.alpha[0,i,1]);



plt.figure()
plt.plot(ff,np.abs(S11))
#plt.plot(ff,np.abs(S12));
plt.show()


eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
#plt.axis('off')
plt.show()

ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()
