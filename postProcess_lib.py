import numpy as np
import math
from array import array
import sys
import matplotlib.pyplot as plt
import os

def readnek( fname ):

	#--------------------------------------------
	# Open the file
	#--------------------------------------------
    arr = array('L')
    with open(fname,'rb') as f:
	
	#--------------------------------------------
        # Read header of tile
    	#--------------------------------------------
    	header = f.readline(132)
 	etag = f.read(4)
 	# Precision of simulation. 
    	wdsz = float(header[5])
    
    	# Element sizes
    	lr1 = [float(filter(None,header[7:9])),float(filter(None,header[10:13])),float(filter(None,header[13:16]))]

    	# Compute total number of points per element
    	npel = int(reduce(lambda x, y: x*y, lr1))

    	# Compute number of active dimensions
    	if (lr1[2] > 1):
    	    add_dims = 1
    	else:
    	    add_dims = 0
    	ndim = 2 + add_dims
    
    	# Number of elements
    	nel = float(filter(None,header[16:26]))
    	
    	# Number of elements in the file
    	nelf = int(filter(None,header[27:37]))
    	
    	# Time
    	time = float(filter(None,header[38:58]))
    	
    	# Iteration number
    	istep = float(filter(None,header[59:68]))
    
    	# Get fields [XUPT]
    	fields = header[83:]
    	var = np.zeros(5)
    	if (int('X' in fields) == 1):
    	    var[0] = ndim
    	if (int('U' in fields) == 1):
    	    var[1] = ndim
    	if (int('P' in fields) == 1):
    	    var[2] = 1
    	if (int('T' in fields) == 1):
    	    var[3] = 1

    	nfields = int(reduce(lambda x, y: x+y, var))
    
    	# Read element map
	#map_proxy = f.readlines()[1:]
    	elmap = np.fromfile(f,dtype='int32',count=nelf)  	
	#--------------------------------------------
        # Read data
        #--------------------------------------------
	data = np.zeros((nelf,npel,nfields))
	for ivar in range(1,6):
	    if (ivar == 1):
		idim0 = 0
	    else: 
		idim0 = reduce(lambda x, y: x+y, var[0:ivar-1])
	    for iel in elmap:
		iter_range = [x+idim0 for x in range(1,int(var[ivar-1])+1)]
		iter_range = [int(x) for x in iter_range]
		for idim in iter_range:
		    data[iel-1,:,idim-1] = np.fromfile(f,dtype='float32',count=npel) 	

    return [data,time,istep]

def reshapenek( data, nelx, nely ):
	nel,N2,nfld = data.shape
	N = math.sqrt(N2)
	if (nel != nelx*nely):
	    print 'Error: nel != nelx*nely.'
	    sys.exit()
	
	#--------------------------------------------
        # Reshape data
        #--------------------------------------------

	mesh = np.zeros((int((N-1)*nelx+1),int((N-1)*nely+1),nfld))

	for ifld in range(0,nfld):
	    # Check this mesh isn't transposed!!!
#	    mesh = np.zeros((int((N-1)*nelx+1),int((N-1)*nely+1)))
	    
	    for iel in range(0,nel):
		ielx = math.floor(iel/nely) + 1
		iely = (iel % nely) + 1
		
		ii = [x+(N-1)*(ielx-1) for x in range(0,int(N))]
		ii = [int(x) for x in ii]
		jj = [x+(N-1)*(iely-1) for x in range(0,int(N))]
		jj = [int(x) for x in jj]

		mesh[ii[0]:(ii[7]+1),jj[0]:(jj[7]+1),ifld] = np.reshape(data[iel,:,ifld], (8,8))
	
	return [ mesh ]


def geometricRatio( r, n, Sn ):
	# Calculating the axis using a geometric ratio.  The variable r is the geometric ratio to be used, n is the number of elements and Sn is the length of the axis, clusterEdge chooses which side to cluster the gridpoints on.

	# Compute first step size.
	a = (r - 1)*Sn/(r**n - 1)
	
	x = None
	geosum = 0
	geosum_last = 0
	
	for i in range(0,n):

	    geosum = geosum + a*r**(i)

	    # Computing spacing for each element
	    small_vector = np.linspace(geosum_last,geosum,8)

	    # Concantenating vectors.
	    if (i < n-1):
		small_vector = small_vector[0:7] 
	    if (i == 0):
		x = small_vector
	    else:
		x = np.concatenate([x, small_vector])
	    geosum_last = geosum

	return [ x ]


def myPcolour(x,y,data,x_label,y_label,x_range,y_range,filename,name,file_counter,**kwargs):
	# Plots figure easily, without having to repeat everything multiple times.
	
        plt.figure(figsize=(50, 50)) # Increases resolution.
        plt.xlabel(x_label,fontsize=80)
        plt.ylabel(y_label,fontsize=80)
        plt.xticks(x_range, fontsize = 60)
        plt.yticks(y_range, fontsize = 60)
        plt.pcolor(x,y,data,**kwargs)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize = 60)  # vertically oriented colorbar

        plt.savefig(''.join([filename,name,repr(file_counter).zfill(5),'.png']))

	plt.close('all')

	return

def particlePcolour(x,y,data,x_label,y_label,x_range,y_range,filename,name,file_counter,x_ppos,y_ppos,**kwargs):
        # Plots figure easily, without having to repeat everything multiple times.

        plt.figure(figsize=(50, 50)) # Increases resolution.
        plt.xlabel(x_label,fontsize=80)
        plt.ylabel(y_label,fontsize=80)
        plt.xticks(x_range, fontsize = 60)
        plt.yticks(y_range, fontsize = 60)
        plt.pcolor(x,y,data,**kwargs)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize = 60)  # vertically oriented colorbar

	plt.scatter(x_ppos,y_ppos,marker='.',color='black')

        plt.savefig(''.join([filename,name,repr(file_counter).zfill(5),'.png']))

        plt.close('all')

        return


def PseudoColourPlotting( filename, jump, total_timesteps, numPlots, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType, particles ):
	# Plots data from a Nek5000 run.  Inputs are as follows:
	# filename: name that comes before the 0.f##### in the output files from Nek5000.
	# jump: number of 0.f##### files to skip between each plot.
	# total_timesteps: number of the last 0.f##### file to consider.
	# numPlots: number of plots to produce (1 - temperature only, 2 - temperature and vertical velocity, 3 - temperature, vertical velocity, and magnitude of velocity).
	# elements_x: number of elements in the x-direction.
	# elements_y: number of elements in the y -direction.
	# gridpoints_x: number of gridpoints in the x-direction.
	# gridpoints_y: number of gridpoints in the y-direction.
	# x_cluster: geometric ratio used to cluster gridpoints in the x-direction.
	# y_cluster: geometric ratio used to cluster gridpoints in the y-direction.
	# gridType: 0 - half domain (i.e. x goes from 0-50 while y goes from 0-100 with a half-plume), 1 - full domain (i.e. domain is square).
	# particles: switch to plot particle paths if particles are included in the simulation.
	
	total_files = total_timesteps/jump;

	file_counter = 1

	range_vals = np.array(range(1,total_files))*jump
	for k in range_vals:

	    # Outputs counter to terminal.
            files_remaining = total_files - k/jump
            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

	    # Reads data files.
	    data,time,istep = readnek(''.join([filename,'0.f',repr(k+1).zfill(5)]))
	    # Reshapes data onto uniform grid.
	    [ mesh ] = reshapenek(data, elements_y, elements_x)
	    # Consider only the necessary number of plots.
	    if (numPlots == 1):
		temperature = mesh[:,:,3]
	    elif (numPlots == 2):
		temperature = mesh[:,:,3]
		verVel = mesh[:,:,1]
	    elif (numPlots == 3):
		verVel = mesh[:,:,1]
		horVel = mesh[:,:,0]
		temperature = mesh[:,:,3]
		magVel = np.sqrt(np.square(verVel) + np.square(horVel))
	    # Defines size of grid.
	    if( gridType == 0 ):
		[ x ] = geometricRatio(1/x_cluster,elements_x,gridpoints_x)
	    else:
	    	[ x1 ] = geometricRatio(1/x_cluster,elements_x/2,gridpoints_x/2)
	    	[ x2 ] = geometricRatio(x_cluster,elements_x/2,gridpoints_x/2)
	    	x = np.concatenate([ x1[:-1], [ x+50 for x in x2 ] ])

	    [ y ] = geometricRatio(y_cluster,elements_y,gridpoints_y)
	    
	    # Reading in particle data.
	    if (particles == 1):
		npart = (k+1)*4
                pname = ''.join(['part',repr(npart+1).zfill(5),'.3D'])
                text_file = open(pname,"r")

                lines = text_file.read().strip()
                lines = lines.splitlines()
                lines = np.asarray([ map(float, line.split()) for line in lines ])
                x_pos = lines[:,0]
                y_pos = lines[:,1]

	    for plotNum in range(0,numPlots):
		if (plotNum == 0):
		    dataPlot = temperature
		    c_min = 0
		    c_max = 2
		    name = '_temp_'
		elif (plotNum == 1):
		    dataPlot = verVel
                    c_min = -1
                    c_max = 5
		    name = '_verVel_'
	    	elif (plotNum == 2):
                    dataPlot = magVel
                    c_min = 0
                    c_max = 5
		    name = '_magVel_'

		if (particles == 0):
		    # Plots reshaped data
		    myPcolour(x,y,dataPlot,'Horizontal position','Vertical position',range(0,101,10),range(0,101,10),filename,name,file_counter,vmin=c_min,vmax=c_max)
		elif (particles == 1):
		    particlePcolour(x,y,dataPlot,'Horizontal position','Vertical position',range(0,101,10),range(0,101,10),filename,name,file_counter,x_pos,y_pos,vmin=c_min,vmax=c_max)



	    file_counter = file_counter + 1

	return

def myPlot(x,y,x_label,y_label,xlimits,ylimits,filename,name,file_counter,**kwargs):
        # Plots figure easily, without having to repeat everything multiple times.

        plt.figure(figsize=(30, 50)) # Increases resolution.
        plt.xlabel(x_label,fontsize=80)
        plt.ylabel(y_label,fontsize=80)
	plt.xlim(xlimits)
	plt.ylim(ylimits)
	plt.tick_params(axis='both', which='major', labelsize=60)
	plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(x,y,**kwargs)

        plt.savefig(''.join([filename,name,repr(file_counter).zfill(5),'.png']))

        plt.close('all')

        return


def integrateDomain( filename, jump, total_timesteps, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType ):
	# Plots line data from a Nek5000 run.  Inputs are as follows:
        # filename: name that comes before the 0.f##### in the output files from Nek5000.
        # jump: number of 0.f##### files to skip between each plot.
        # total_timesteps: number of the last 0.f##### file to consider.
        # elements_x: number of elements in the x-direction.
        # elements_y: number of elements in the y -direction.
        # gridpoints_x: number of gridpoints in the x-direction.
        # gridpoints_y: number of gridpoints in the y-direction.
        # x_cluster: geometric ratio used to cluster gridpoints in the x-direction.
        # y_cluster: geometric ratio used to cluster gridpoints in the y-direction.
        # gridType: 0 - half domain (i.e. x goes from 0-50 while y goes from 0-100 with a half-plume), 1 - full domain (i.e. domain is square).

	total_files = total_timesteps/jump;
	ambientTemp = 273.15;
	g = 9.81

	file_counter = 1

        range_vals = np.array(range(0,total_files))*jump
        for k in range_vals:

	    # Outputs counter to terminal.
            files_remaining = total_files - k/jump
            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            data,time,istep = readnek(''.join([filename,'0.f',repr(k+1).zfill(5)]))
            # Reshapes data onto uniform grid.
            [ mesh ] = reshapenek(data, elements_y, elements_x)
            verVel = mesh[:,:,3]
	    horVel = mesh[:,:,2]
	    #magVel = np.sqrt(np.square(verVel) + np.square(horVel))
            temperature = mesh[:,:,5]
	    temperature = temperature + 273.15
            # Defines size of grid.
	    y_length,x_length = verVel.shape
	    
	    # Defines size of grid.
            if( gridType == 0 ):
                [ x ] = geometricRatio(1/x_cluster,elements_x,gridpoints_x)
            else:
                [ x1 ] = geometricRatio(1/x_cluster,elements_x/2,gridpoints_x/2)
                [ x2 ] = geometricRatio(x_cluster,elements_x/2,gridpoints_x/2)
                x = np.concatenate([ x1[:-1], [ x+50 for x in x2 ] ])

            [ y ] = geometricRatio(y_cluster,elements_y,gridpoints_y)



	    buoyancy_flux = np.zeros((y_length))
	    velocity_int = np.zeros((y_length))
	    for height in range(0,y_length):
		local_buoy = np.zeros((x_length))
		for i in range(0,x_length):
		    # Calculate w*g' at each point on the horizontal
		    #if (temperature[height,i] == 0):
                    #    reduced_gravity = 9.81
		    #else:
		    reduced_gravity = g*(temperature[height,i] - ambientTemp)/temperature[height,i]
		    local_buoy[i] = reduced_gravity*verVel[height,i]
		    #if (int(math.isnan(local_buoy[k])) == 1):
		#	local_buoy[k] = 0
		for i in range(1,x_length-1):
		    dx = x[i] - x[i-1]
		    buoyancy_flux[height] = buoyancy_flux[height] + (dx/6)*(local_buoy[i-1] + 4*local_buoy[i] + local_buoy[i+1])
		    velocity_int[height] = velocity_int[height] + (dx/6)*(verVel[height,i-1] + 4*verVel[height,i] + verVel[height,i+1])

	    # Plots integral data at each time-step.

	    myPlot(buoyancy_flux,y,'Buoyancy','Vertical position',[0,0.08],[0,100],filename,'_buoy_',file_counter)
	    myPlot(velocity_int,y,'Vertical Velocity','Vertical position',[0,12],[0,100],filename,'_velInt__',file_counter)

            file_counter = file_counter + 1

	return

# The function PseudoColourPlotting takes the following form and plots data from a Nek5000 run.
# PseudoColourPlotting( filename, jump, total_timesteps, numPlots, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType )
# Inputs are as follows:
        # filename: name that comes before the 0.f##### in the output files from Nek5000.
        # jump: number of 0.f##### files to skip between each plot.
        # total_timesteps: number of 0.f##### files to consider (not number of last file).
        # numPlots: number of plots to produce (1 - temperature only, 2 - temperature and vertical velocity, 3 - temperature, vertical velocity, and magnitude of velocity).
        # elements_x: number of elements in the x-direction.
        # elements_y: number of elements in the y -direction.
        # gridpoints_x: number of gridpoints in the x-direction.
        # gridpoints_y: number of gridpoints in the y-direction.
        # x_cluster: geometric ratio used to cluster gridpoints in the x-direction.
        # y_cluster: geometric ratio used to cluster gridpoints in the y-direction.
        # gridType: 0 - half domain (i.e. x goes from 0-50 while y goes from 0-100 with a half-plume), 1 - full domain (i.e. domain is square).

#PseudoColourPlotting( 'plume_v2_largeIni', 1, 1, 3, 40, 50, 100, 100, 1.25, 1.1, 1 )


# The function integrateDomain takes the following form and plots line data from a Nek5000 run.
# integrateDomain( filename, jump, total_timesteps, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType )
# Inputs are as follows:
        # filename: name that comes before the 0.f##### in the output files from Nek5000.
        # jump: number of 0.f##### files to skip between each plot.
        # total_timesteps: number of 0.f##### files to consider (not number of last file).
        # elements_x: number of elements in the x-direction.
        # elements_y: number of elements in the y -direction.
        # gridpoints_x: number of gridpoints in the x-direction.
        # gridpoints_y: number of gridpoints in the y-direction.
        # x_cluster: geometric ratio used to cluster gridpoints in the x-direction.
        # y_cluster: geometric ratio used to cluster gridpoints in the y-direction.
        # gridType: 0 - half domain (i.e. x goes from 0-50 while y goes from 0-100 with a half-plume), 1 - full domain (i.e. domain is square).

#integrateDomain( 'plume_v2_largeIni', 10, 2400, 40, 80, 50, 100 )

# Combine images into video using the script make_video_from_image.sh
#os.system('./make_video_from_image.sh')


