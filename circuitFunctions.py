import meep as mp

def createCircuitGeometry(a, resolution, mwavelength, additional_buffer):

	#need to have the PML thickness equal to 1/2 the maximum wavelength
	#to allow for wide bandwidth, pick this as just the maximum wavelength
	pmlThickness = mwavelength/(2*a);
	
	if pmlThickness < 1:
		pmlThickness = 1
		
	print("PML Thickness: ", pmlThickness)
	buffer = 5/resolution + additional_buffer;
	
	pml_layers = [mp.PML(pmlThickness)]
	
	#setup the board 
	bwidth = 5*.0254/a;
	bheight = 4*.0254/a;
	bthick = mp.inf;
	bdielectric = 4.4;
	
	board = mp.Block(mp.Vector3(bwidth,bheight,bthick),
                     center=mp.Vector3(0,0,0),
                     material=mp.Medium(epsilon=bdielectric));

	#setup the trace
	twidth = 3.5*.0254/a;
	theight = 0.2*.0254/a;
	tthick = mp.inf;

	trace = mp.Block(mp.Vector3(twidth,theight,tthick),
                   center=mp.Vector3(0,0,0),
                   material=mp.metal);  
                   
	geometry = [board,trace];
	cell = mp.Vector3(bwidth+2*(buffer+pmlThickness),bheight+2*(buffer+pmlThickness),0);
	
	#setup the near field box
	nfXPos = 0.5*(bwidth+buffer/2);
	nfYPos =0.5*(bheight+buffer/2);
	nfXWidth = mp.Vector3(x=bwidth);
	nfYWidth = mp.Vector3(y=bheight);
	
	#initialize the port location
	pWidth = 2/resolution
	p1Loc = mp.Vector3(-twidth/2,0,0);
	
	nPointsPerTraceThickness = resolution/theight;
	
	print("Points Per Trace Thickness: ", nPointsPerTraceThickness)
	
	return geometry, cell, pml_layers, p1Loc , nfXPos, nfYPos, nfXWidth, nfYWidth, board, buffer,pmlThickness           