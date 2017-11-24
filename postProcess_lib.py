import numpy as np
import math
from array import array
import sys
import matplotlib.pyplot as plt
# The following allows the script to work over SSH.  Check it still works when not SSH tunnelling.
plt.switch_backend('agg')
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


def myPcolour(x,y,data,time,x_label,y_label,x_range,y_range,filename,name,file_counter,**kwargs):
	# Plots figure easily, without having to repeat everything multiple times.
	
        plt.figure(figsize=(25, 25)) # Increases resolution.
	plt.title('time = %d'%(time),fontsize=40)
        plt.xlabel(x_label,fontsize=40)
        plt.ylabel(y_label,fontsize=40)
        plt.xticks(x_range, fontsize = 30)
        plt.yticks(y_range, fontsize = 30)
        plt.pcolormesh(x,y,data,**kwargs)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize = 30)  # vertically oriented colorbar

        plt.savefig(''.join([filename,name,repr(file_counter).zfill(5),'.png']))

	plt.close('all')

	return

def particlePcolour(x,y,data,time,x_label,y_label,x_range,y_range,filename,name,file_counter,x_ppos,y_ppos,**kwargs):
        # Plots figure easily, without having to repeat everything multiple times.

	particlesOnly = 0
	
	if(particlesOnly < 1):
	    plt.figure(figsize=(25, 25)) # Increases resolution.
	    plt.title('time = %s'%round(time,2),fontsize=40)
            plt.xlabel(x_label,fontsize=40)
            plt.ylabel(y_label,fontsize=40)
            plt.xticks(x_range, fontsize = 30)
            plt.yticks(y_range, fontsize = 30)
            plt.pcolormesh(x,y,data,**kwargs)
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize = 30)  # vertically oriented colorbar

	plt.scatter(x_ppos,y_ppos,marker='.',color='black',s=50)

	#if(particlesOnly > 0):
	plt.axis([0,50,0,100])

        plt.savefig(''.join([filename,name,repr(file_counter).zfill(5),'_particle.png']))

        plt.close('all')

        return


def PseudoColourPlotting( filename, start_file, jump, total_timesteps, numPlots, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType, particles ):
	# Plots data from a Nek5000 run.  Inputs are as follows:
	# filename: name that comes before the 0.f##### in the output files from Nek5000.
	# start_file: file number to start at (usually leave at 1).
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

	#file_counter = 1
	
	if (start_file == 1):
	    range_vals = [x - (jump - 1) for x in np.array(range(1,total_files))*jump]
	else:
	    range_vals = np.array(range(int(math.floor(start_file/jump)),total_files))*jump
	    print "Make sure calculation of file numbers is correct"
	    print range_vals

	for k in range_vals:

	    # Outputs counter to terminal.
	    if (start_file == 1):
		files_remaining = total_files - k/jump
	    else:
            	files_remaining = total_files - k/jump - start_file/jump

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

	    # Reads data files.
	    data,time,istep = readnek(''.join([filename,'0.f',repr(k).zfill(5)]))
	    # Reshapes data onto uniform grid.
	    [ mesh ] = reshapenek(data, elements_y, elements_x)
	    # Consider only the necessary number of plots.
	    if (particles == 1):
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
	    else:
		if (numPlots == 1):
                    temperature = mesh[:,:,5]
                elif (numPlots == 2):
                    temperature = mesh[:,:,5]
                    verVel = mesh[:,:,3]
                elif (numPlots == 3):
                    verVel = mesh[:,:,3]
                    horVel = mesh[:,:,2]
                    temperature = mesh[:,:,5]
                    magVel = np.sqrt(np.square(verVel) + np.square(horVel))

	    # Defines size of grid.
	    if( gridType == 0 ):
		[ x ] = geometricRatio(x_cluster,elements_x,gridpoints_x)
	    else:
	    	[ x1 ] = geometricRatio(x_cluster,elements_x/2,gridpoints_x/2)
	    	[ x2 ] = geometricRatio(1/x_cluster,elements_x/2,gridpoints_x/2)
	    	x = np.concatenate([ x1[:-1], [ x+50 for x in x2 ] ])

	    [ y ] = geometricRatio(y_cluster,elements_y,gridpoints_y)
	    
	    # Reading in particle data.
	    if (particles == 1):
		npart = (k)*4
                pname = ''.join(['part',repr(npart).zfill(5),'.3D'])
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
		    c_max = 0.3
		    name = '_temp_'
		elif (plotNum == 1):
		    dataPlot = verVel
                    c_min = 0
                    c_max = 10
		    name = '_verVel_'
	    	elif (plotNum == 2):
                    dataPlot = magVel
                    c_min = 0
                    c_max = 10
		    name = '_magVel_'

		if (particles == 0):
		    # Plots reshaped data
		    myPcolour(x,y,dataPlot,time,'Horizontal position','Vertical position',range(0,101,10),range(0,101,10),filename,name,file_counter,vmin=c_min,vmax=c_max)
		elif (particles == 1):
		    particlePcolour(x,y,dataPlot,time,'Horizontal position','Vertical position',range(0,101,10),range(0,101,10),filename,name,k/jump-1,x_pos,y_pos,vmin=c_min,vmax=c_max)

	    #file_counter = file_counter + 1

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

def meshInd( filename, jump, total_timesteps, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType ):
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

	total_files = total_timesteps/jump
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
	    mean_u = np.mean(verVel)
	    mean_v = np.mean(horVel)

	    u_prime = [(x-mean_u)**2 for x in verVel]
	    v_prime = [(x-mean_v)**2 for x in horVel]

	    k = 0.5*(np.mean(u_prime) + np.mean(v_prime))
	    print k

	return


def particleAdvect( filename, jump, total_timesteps, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType ):
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

	# Initial values of particle parameters.
        initial_particle_x = 49.8
        initial_particle_y = 0.2
        particle_x_velocity = 0.0
        particle_y_velocity = 0.0

	# Calulation of settling velocity of particle.
	fluid_viscosity = 8.76e-6
	fluid_density = 1000 
	particle_density = 2000
	particle_diameter = 0.0001
	particle_settlingVel = 0 #9.81*particle_diameter**2*(particle_density - fluid_density)/(18*fluid_viscosity)
	print "Particle Settling Velocity: {:2d}".format(particle_settlingVel)

        # Initialisation of loop.
        file_counter = 1
        time_old = 0
        range_vals = np.array(range(1,total_files))*jump
	particle_position_x = initial_particle_x
	particle_position_y = initial_particle_y
	particle_position_vector = np.zeros((total_files,2))
	particle_position_vector[0,:] = [particle_position_x,particle_position_y]

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
	  
	    # Compute timestep.
	    dt = time - time_old
	    time_old = time

            # Defines size of grid.
            if( gridType == 0 ):
                [ x ] = geometricRatio(x_cluster,elements_x,gridpoints_x)
            else:
                [ x1 ] = geometricRatio(x_cluster,elements_x/2,gridpoints_x/2)
                [ x2 ] = geometricRatio(1/x_cluster,elements_x/2,gridpoints_x/2)
                x = np.concatenate([ x1[:-1], [ x+gridpoints_x/2 for x in x2 ] ])

	    [ y ] = geometricRatio(y_cluster,elements_y,gridpoints_y)

	    #plt.pcolor(x,y,verVel)
	    #plt.savefig(''.join([filename,'_particle2','.png']))

	    # Computes gridbox in which the particle lies.
	    x_range = np.array(range(0,len(x)))
            y_range = np.array(range(0,len(y)))
	    breaker_x = 0
	    breaker_y = 0
	    for x_pos in x_range:
		if(breaker_x < 1):
                    if(x[x_pos] > particle_position_x):
                        i = x_pos
                        breaker_x = 1
	    for y_pos in y_range:
		if(breaker_y < 1):
                    if(y[y_pos] > particle_position_y):
                        j = len(y) - 1 - y_pos
			breaker_y = 1

	    # Computes weights, deciding 'how much' of the velocity from each node surrounding the particle should be transferred to it.
	    w1 = 0.25
	    w2 = 0.25
	    w3 = 0.25
	    w4 = 0.25

	    # Computes velocity of the particle.
	    particle_x_velocity = particle_x_velocity + w1*horVel[j,i] + w2*horVel[j-1,i] + w3*horVel[j,i-1] + w4*horVel[j-1,i-1]
	    particle_y_velocity = particle_y_velocity + w1*verVel[j,i] + w2*verVel[j-1,i] + w3*verVel[j,i-1] + w4*verVel[j-1,i-1]

	    # Advects the particle.
	    particle_position_x = particle_position_x + particle_x_velocity*dt
            particle_position_y = particle_position_y + particle_y_velocity*dt - particle_settlingVel*dt

	    if(particle_position_x > gridpoints_x):
		particle_position_x = gridpoints_x
		particle_x_velocity = 0
	    if(particle_position_x < 0):
                particle_position_x = 0
                particle_x_velocity = 0
            if(particle_position_y > gridpoints_y):
                particle_position_y = gridpoints_y
                particle_y_velocity = 0
            if(particle_position_y < 0):
                particle_position_y = gridpoints_y
                particle_y_velocity = 0

	    particle_position_vector[file_counter,:] = [particle_position_x,particle_position_y]

	    file_counter = file_counter + 1

	print particle_position_vector

	for plot_count in range(0,total_files,100):
	    plt.scatter(particle_position_vector[plot_count,0],particle_position_vector[plot_count,1],marker='.',color='black',s=0.5)
	    axes = plt.gca()
	    axes.set_xlim([0,gridpoints_x])
	    axes.set_ylim([0,gridpoints_y])
	    print "done"
	plt.savefig(''.join([filename,'_particle','.png']))

	return




def meshPlot( elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType ):

	# Defines size of grid.
        if( gridType == 0 ):
            [ x ] = geometricRatio(x_cluster,elements_x,gridpoints_x)
        else:
            [ x1 ] = geometricRatio(x_cluster,elements_x/2,gridpoints_x/2)
            [ x2 ] = geometricRatio(1/x_cluster,elements_x/2,gridpoints_x/2)
            x = np.concatenate([ x1[:-1], [ x+gridpoints_x/2 for x in x2 ] ])

        [ y ] = geometricRatio(y_cluster,elements_y,gridpoints_y)

	for i in range(0,len(x)):
            xplot = np.zeros(len(y))
            xplot = [q + x[i] for q in xplot]
            plt.plot(xplot,y,color='black',linewidth=0.5)
        for j in range(0,len(y)):
            yplot = np.zeros(len(x))
            yplot = [p + y[j] for p in yplot]
            plt.plot(x,yplot,color='black',linewidth=0.5)

        plt.savefig('mesh.jpg')

	return



