import meep as mp

def createCircuitGeometry(a, resolution):

	pmlThickness = 1;
	buffer = 4;
	
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
	theight = 0.3*.0254/a;
	tthick = mp.inf;

	trace = mp.Block(mp.Vector3(twidth,theight,tthick),
                   center=mp.Vector3(0,0,0),
                   material=mp.metal);  
                   
	geometry = [board,trace];
	cell = mp.Vector3(bwidth+buffer,bheight+buffer,0);
	
	#setup the near field box
	nfXPos = 0.5*(bwidth+buffer/2);
	nfYPos =0.5*(bheight+buffer/2);
	nfXWidth = mp.Vector3(x=bwidth);
	nfYWidth = mp.Vector3(y=bheight);
	
	#initialize the port location
	pWidth = 1/resolution
	p1Loc = mp.Vector3(-twidth/2 ,0,0);
	
	nPointsPerTraceThickness = resolution/theight;
	
	print("Points Per Trace Thickness: ", nPointsPerTraceThickness)
	
	return geometry, cell, pml_layers, p1Loc , nfXPos, nfYPos, nfXWidth, nfYWidth           