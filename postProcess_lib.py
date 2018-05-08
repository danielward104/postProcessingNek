import numpy as np
import math
from array import array
import sys
import matplotlib.pyplot as plt
# The following allows the script to work over SSH.  Check it still works when not SSH tunnelling.
plt.switch_backend('agg')
import os
from tempfile import TemporaryFile
from colorsys import hsv_to_rgb

# Import user-defined modules.
import readingNek as rn
import mappings as mp
import plottingTools as pt
import generalTools as tools

def PseudoColourPlotting( filename, order, dimension, start_file, jump, final_timestep, numPlots, elements_x, elements_y, elements_z, x_start, x_end, y_start, y_end, z_slice, x_cluster, y_cluster,  particles ):
	# Plots data from a Nek5000 run.  Inputs are as follows:
	# filename: name that comes before the 0.f##### in the output files from Nek5000.
	# Order of the polynomial used for the spectral element method.
	# start_file: file number to start at (usually leave at 1).
	# jump: number of 0.f##### files to skip between each plot.
	# final_timestep: number of the last 0.f##### file to consider.
	# numPlots: number of plots to produce (1 - temperature only, 2 - temperature and vertical velocity, 3 - temperature, vertical velocity, and magnitude of velocity).
	# elements_x: number of elements in the x-direction.
	# elements_y: number of elements in the y -direction.
	# gridpoints_x: number of gridpoints in the x-direction.
	# gridpoints_y: number of gridpoints in the y-direction.
	# x_cluster: geometric ratio used to cluster gridpoints in the x-direction.
	# y_cluster: geometric ratio used to cluster gridpoints in the y-direction.
	# particles: switch to plot particle paths if particles are included in the simulation.

	quiver = 0	

	final_file = final_timestep/jump;

	#file_counter = 1
	
	if (start_file == 1):
	    range_vals = [x - (jump - 1) for x in np.array(range(1,final_file+1))*jump]

	else:
	    range_vals = np.array(range(int(math.floor(float(start_file)/jump)),final_file+1))*jump
	    print "Make sure calculation of file numbers is correct"
	    print range_vals

	for k in range_vals:

	    # Outputs counter to terminal.
	    if (start_file == 1):
		files_remaining = final_file + 1 - k/jump
	    else:
            	files_remaining = final_file + 1 - (k-range_vals[0])/jump - start_file/jump

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

	    # Reads data files.
	    data,time,istep,header,elmap = rn.readnek(''.join([filename,'0.f',repr(k).zfill(5)]))
	    # Reshapes data onto uniform grid.

	    if (dimension == 2):

	        [ mesh ] = rn.reshapenek2D(data, elements_y, elements_x)

	    elif (dimension == 3):

	        [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)
		mesh = mesh[:,:,4,:]

	    # Consider only the necessary number of plots.
	    # Different numbering system for particles and no-particles.
	
	    if (particles == 1):
	    	if (numPlots == 1):
		    temperature = mesh[:,:,3]

		    if (quiver == 1):
			verVel = mesh[:,:,1]
                    	horVel = mesh[:,:,0]
	
	    	elif (numPlots == 2):
		    verVel = mesh[:,:,1]
		    horVel = mesh[:,:,0]
		    magVel = np.sqrt(np.square(verVel) + np.square(horVel))
		    temperature = mesh[:,:,3]

	    elif (particles == 0):
	        if (numPlots == 1):
                   temperature = np.transpose(mesh[:,:,7])

		    
		   if (quiver == 1):
                        verVel = mesh[:,:,3]
                        horVel = mesh[:,:,2]

                elif (numPlots == 2):
                    verVel = mesh[:,:,3]
                    horVel = mesh[:,:,2]
                    temperature = mesh[:,:,5]
                    magVel = np.sqrt(np.square(verVel) + np.square(horVel))

	    # Defines size of grid.
	    x = mp.mesh_generation(x_cluster,elements_x,x_start,x_end,order,2,'in')
	    y = mp.mesh_generation(y_cluster,elements_y,y_start,y_end,order,2,'in')

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
		    # Plots the square root of temperature to improve the outline of the plume.
		    dataPlot = temperature #[ abs(t)**0.5 for t in temperature ]
		    #c_min = -5
		    #c_max = 5
		    name = 'Temperature'

		    if (quiver == 1):
			quiverPlotx = horVel
			quiverPloty = verVel 

		elif (plotNum == 1):
		    dataPlot = verVel
                    c_min = -1
                    c_max = 1
		    name = 'y-velocity'
	    	elif (plotNum == 2):
                    dataPlot = magVel
                    c_min = 0
                    c_max = 1.0
		    name = 'Velocity-magnitude'
                elif (plotNum == 3):
                    dataPlot = horVel
                    c_min = -0.2
                    c_max = 0.2
                    name = 'x-velocity'
		elif (plotNum == 4):
		    dataPlot = pressure
                    c_min = -0.1
                    c_max = 0.1
                    name = 'Pressure'

		zoom = 0

		if (zoom == 1):
		    x_start = int(2*(len(x)+1)/6)
		    x_end = int(4*(len(x)+1)/6)
		    y_start = int(4*(len(y)+1)/6)
		    y_end = int(6*(len(y)+1)/6)

		    dataPlot = np.array(dataPlot)
		    dataPlot = dataPlot[y_start:y_end,x_start:x_end]

		    x = x[x_start:x_end]
		    y = y[y_start:y_end]
		    name = 'Temperature_close'

		    if (quiver == 1):
			quiverPlotx = quiverPlotx[y_start:y_end,x_start:x_end]
		    	quiverPloty = quiverPloty[y_start:y_end,x_start:x_end]

		if (quiver == 1):
		    pt.myPcolourQuiver(np.transpose(x),y,dataPlot,quiverPlotx,quiverPloty,time,'Horizontal position','Vertical position',filename,name,k/jump,vmin=c_min,vmax=c_max,cmap='RdBu_r')
		else:
		    if (particles == 0):
		    	# Plots reshaped data
#		    	pt.myPcolour(np.transpose(x),y,dataPlot,time,'Horizontal position','Vertical position',filename,name,k/jump,vmin=c_min,vmax=c_max,cmap='RdBu_r')
			pt.myPcolour(np.transpose(x),y,dataPlot,time,'Horizontal position','Vertical position',filename,name,k/jump,cmap='RdBu_r')
		    elif (particles == 1):
		    	pt.myPcolour(y,np.transpose(x),dataPlot,time,'Horizontal position','Vertical position',filename,name,k/jump,vmin=c_min,vmax=c_max,cmap='RdBu_r')

#		    pt.particlePcolour(np.transpose(x),y,dataPlot,time,'Horizontal position','Vertical position',range(0,51,5),range(0,101,10),filename,name,k/jump,x_pos,y_pos,vmin=c_min,vmax=c_max,cmap='RdBu_r')

	return


def trapezium( points_x, points_y, data ):

	# Computes the approximate two-dimensional integral of the function represented by 'data'.

	sum_total = 0
	x_tot = np.shape(points_x)
	x_tot = x_tot[0] - 1
	y_tot = np.shape(points_y)
	y_tot = y_tot[0] - 1

	for x in range(0,x_tot):
	    for y in range(0,y_tot):
		sum1 = data[y,x] + data[y+1,x] + data[y,x+1] + data[y+1,x+1]
		dx = points_x[x+1] - points_x[x]
		dy = points_y[y+1] - points_y[y]
		trap = sum1*dx*dy/4
		sum_total = sum_total + trap

	return sum_total


def integrateDomain( filename, order, jump, total_timesteps, elements_x, elements_y, x_start, x_end, y_start, y_end, x_cluster, y_cluster, gridType ):

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
	ambientTemp = 273.15
	ambientDensity = 1027
	g = 9.81
	expCoeff = 0.0002

	file_counter = 1

        range_vals = np.array(range(0,total_files))*jump
	buoyancy_flux = np.zeros(total_files+1)
	time_vec = np.zeros(total_files+1)
	all_energies = np.zeros(total_files+1)
	all_buoyancy = np.zeros(total_files+1)
	all_avgVels = np.zeros(total_files+1)

        for k in range_vals:

	    # Outputs counter to terminal.
            files_remaining = total_files - k/jump
            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            data,time,istep = rn.readnek(''.join([filename,'0.f',repr(k+1).zfill(5)]))
            # Reshapes data onto uniform grid.
            [ mesh ] = rn.reshapenek(data, elements_y, elements_x)
            verVel = mesh[:,:,3]
	    horVel = mesh[:,:,2]
	    magVel = np.sqrt(np.square(verVel) + np.square(horVel))
            temperature = mesh[:,:,5]
	    temperature = temperature + 273.15
            # Defines size of grid.
	    y_length,x_length = temperature.shape

	    time_vec[file_counter] = time

	    # Defines size of grid.
	    x = mp.mesh_generation(x_cluster,elements_x,x_start,x_end,order,2,'in')
            y = mp.mesh_generation(y_cluster,elements_y,y_start,y_end,order,2,'in')

	    # Computing the integral of the energy and buoyancy.

	    density = [1027*(2+273)/T for T in temperature]
	    density = [ambientDensity*(1-expCoeff*(T - ambientTemp)) for T in temperature]
	    energy = np.square(magVel)
	    energy = 0.5*np.multiply(density,energy)
	    energy_at_time = trapezium(x,y,energy)
	    buoyancy = [g*(p-ambientDensity)/p for p in density]
	    buoyancy_total = trapezium(x,y,np.array(buoyancy))
	    
	    all_buoyancy[file_counter] = buoyancy_total
	    all_energies[file_counter] = energy_at_time

	    # Average velocity.

	    avgVel = np.sum(np.sum(magVel))/(x_length*y_length)
	    avgVel = trapezium(x,y,magVel)
	    all_avgVels[file_counter] = avgVel

	    file_counter = file_counter + 1


        vec1 = np.stack((time_vec, all_energies),axis=1)
        f = open("kinetic_energy.txt","w")
	f.write("# Time, Kinetic Energy`\n")
        np.savetxt(f, vec1)
	f.close()

	vec2 = np.stack((time_vec, all_buoyancy),axis=1)
        f = open("buoyancy.txt","w")
        f.write("# Time, Buoyancy`\n")
        np.savetxt(f, vec2)
        f.close()

	vec3 = np.stack((time_vec, all_avgVels),axis=1)
        f = open("avgVels.txt","w")
        f.write("# Time, Average Velocity`\n")
        np.savetxt(f, vec3)
        f.close()

	return

 
def meshInd():

    	time_4020, ke_4020 = np.loadtxt("./save_v1_40x20/kinetic_energy.txt", unpack = True)
        time_8040, ke_8040 = np.loadtxt("./save_v2_80x40/kinetic_energy.txt", unpack = True)
	time_12060, ke_12060 = np.loadtxt("./save_v3_120x60/kinetic_energy.txt", unpack = True)
	time_16080, ke_16080 = np.loadtxt("./save_v4_160x80/kinetic_energy.txt", unpack = True)

        time_4020, buoy_4020 = np.loadtxt("./save_v1_40x20/buoyancy.txt", unpack = True)
        time_8040, buoy_8040 = np.loadtxt("./save_v2_80x40/buoyancy.txt", unpack = True)
        time_12060, buoy_12060 = np.loadtxt("./save_v3_120x60/buoyancy.txt", unpack = True)
	time_16080, buoy_16080 = np.loadtxt("./save_v4_160x80/buoyancy.txt", unpack = True)

        time_4020, avgVel_4020 = np.loadtxt("./save_v1_40x20/avgVels.txt", unpack = True)
        time_8040, avgVel_8040 = np.loadtxt("./save_v2_80x40/avgVels.txt", unpack = True)
        time_12060, avgVel_12060 = np.loadtxt("./save_v3_120x60/avgVels.txt", unpack = True)
        time_16080, avgVel_16080 = np.loadtxt("./save_v4_160x80/avgVels.txt", unpack = True)

	# Finds time such that all simulations have been run to the same (ish) time.
	min_time =  min(max(time_4020),max(time_8040),max(time_12060),max(time_16080))
	length_4020 = tools.find_nearest(time_4020,min_time)
	length_8040 = tools.find_nearest(time_8040,min_time)
	length_12060 = tools.find_nearest(time_12060,min_time)
	length_16080 = tools.find_nearest(time_16080,min_time)

	# Plots kinetic energy vs. time.
        plt.figure(figsize=(50, 30)) # Increases resolution.
	ax = plt.axes()
	plt.xlabel('Time',fontsize=80)
	plt.ylabel('Kinetic Energy',fontsize=80)
	plt.tick_params(axis='both', which='major', labelsize=60)
	plt.tick_params(axis='both', which='minor', labelsize=60)
	plt.plot(time_4020[2:length_4020],ke_4020[2:length_4020], label="Grid 40x20", linewidth = 5.0)
	plt.plot(time_8040[2:length_8040],ke_8040[2:length_8040], label="Grid 80x40", linewidth = 5.0)
	plt.plot(time_12060[2:length_12060],ke_12060[2:length_12060], label="Grid 120x60", linewidth = 5.0)
        plt.plot(time_16080[2:length_16080],ke_16080[2:length_16080], label="Grid 160x80", linewidth = 5.0)
	ax.yaxis.get_offset_text().set_fontsize(50)
	plt.legend(fontsize=40)
	plt.savefig(''.join(['plume_v3_stratified_keTime.png']))    
	plt.close('all')

        # Plots buoyancy vs. time.
        plt.figure(figsize=(50, 30)) # Increases resolution.
        ax = plt.axes()
        plt.xlabel('Time',fontsize=80)
        plt.ylabel('Buoyancy',fontsize=80)
        plt.tick_params(axis='both', which='major', labelsize=60)
        plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(time_4020[2:length_4020],buoy_4020[2:length_4020], label="Grid 40x20", linewidth = 5.0)
        plt.plot(time_8040[2:length_8040],buoy_8040[2:length_8040], label="Grid 80x40", linewidth = 5.0)
        plt.plot(time_12060[2:length_12060],buoy_12060[2:length_12060], label="Grid 120x60", linewidth = 5.0)
	plt.plot(time_16080[2:length_16080],buoy_16080[2:length_16080], label="Grid 160x80", linewidth = 5.0)
        ax.yaxis.get_offset_text().set_fontsize(50)
        plt.legend(fontsize=40)
        plt.savefig(''.join(['plume_v3_stratified_buoyTime.png']))
        plt.close('all')

	# Plots average velocity.

        plt.figure(figsize=(50, 30)) # Increases resolution.
        ax = plt.axes()
        plt.xlabel('Time',fontsize=80)
        plt.ylabel('Average Velocity',fontsize=80)
        plt.tick_params(axis='both', which='major', labelsize=60)
        plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(time_4020[2:length_4020],avgVel_4020[2:length_4020], label="Grid 40x20", linewidth = 5.0)
        plt.plot(time_8040[2:length_8040],avgVel_8040[2:length_8040], label="Grid 80x40", linewidth = 5.0)
        plt.plot(time_12060[2:length_12060],avgVel_12060[2:length_12060], label="Grid 120x60", linewidth = 5.0)
        plt.plot(time_16080[2:length_16080],avgVel_16080[2:length_16080], label="Grid 160x80", linewidth = 5.0)
        ax.yaxis.get_offset_text().set_fontsize(50)
        plt.legend(fontsize=40)
        plt.savefig(''.join(['plume_v3_stratified_avgVel.png']))
        plt.close('all')


	# Total kinetic energy of the system computed using the trapezium rule.
	tot_ke_4020 = 0
	tot_ke_8040 = 0
	tot_ke_12060 = 0
	for i in range(1,length_4020-1):
	    tot_ke_4020 = tot_ke_4020 + (time_4020[i+2]-time_4020[i+1])*(ke_4020[i+1]+ke_4020[i+2])/2
	for i in range(1,length_8040-1):
            tot_ke_8040 = tot_ke_8040 + (time_8040[i+2]-time_8040[i+1])*(ke_8040[i+1]+ke_8040[i+2])/2
	for i in range(1,length_12060-1):
            tot_ke_12060 = tot_ke_12060 + (time_12060[i+2]-time_12060[i+1])*(ke_12060[i+1]+ke_12060[i+2])/2	

	total_ke = [tot_ke_4020,tot_ke_8040,tot_ke_12060]

	plt.figure(figsize=(50, 30)) # Increases resolution.
        ax = plt.axes()
        plt.xlabel('Number of Elements',fontsize=80)
        plt.ylabel('Total Kinetic Energy',fontsize=80)
        plt.tick_params(axis='both', which='major', labelsize=60)
        plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(elements,total_ke,linewidth = 5.0)
        ax.yaxis.get_offset_text().set_fontsize(50)
        plt.savefig(''.join(['plume_v3_stratified_keElement.png']))

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
            data,time,istep = rn.readnek(''.join([filename,'0.f',repr(k+1).zfill(5)]))
            # Reshapes data onto uniform grid.
            [ mesh ] = rn.reshapenek(data, elements_y, elements_x)
	    verVel = mesh[:,:,1]
            horVel = mesh[:,:,0]
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
                        j = y_pos #len(y) - 1 - y_pos
			breaker_y = 1

	    # Computes weights, deciding 'how much' of the velocity from each node surrounding the particle should be transferred to it.

	    xp = particle_position_x
	    yp = particle_position_y
	    area = (x[i] - x[i-1])*(y[j]-y[j-1])

	    w1 = (xp - x[i-1])*(yp - y[j-1])/area
	    w2 = (xp - x[i-1])*(y[j] - yp)/area
	    w3 = (x[i] - xp)*(yp - y[j-1])/area
	    w4 = (x[i] - xp)*(y[j] - yp)/area

	    j = len(y) - 1 - j

	    # Computes velocity of the particle.
	    particle_x_velocity = w1*horVel[j,i] + w2*horVel[j-1,i] + w3*horVel[j,i-1] + w4*horVel[j-1,i-1]
	    particle_y_velocity = w1*verVel[j,i] + w2*verVel[j-1,i] + w3*verVel[j,i-1] + w4*verVel[j-1,i-1]

	    # Advects the particle.
	    particle_position_x = particle_position_x + particle_x_velocity*dt
            particle_position_y = particle_position_y + particle_y_velocity*dt #- particle_settlingVel*dt

#	    if(particle_position_x > gridpoints_x):
#		particle_position_x = gridpoints_x
#		particle_x_velocity = 0
#	    if(particle_position_x < 0):
#                particle_position_x = 0
#                particle_x_velocity = 0
#            if(particle_position_y > gridpoints_y):
#                particle_position_y = gridpoints_y
#                particle_y_velocity = 0
#            if(particle_position_y < 0):
#                particle_position_y = gridpoints_y
#                particle_y_velocity = 0

	    particle_position_vector[file_counter,:] = [particle_position_x,particle_position_y]

	    file_counter = file_counter + 1

	print particle_position_vector

	for plot_count in range(0,total_files,1):
	    plt.scatter(particle_position_vector[plot_count,0],particle_position_vector[plot_count,1],marker='.',color='black',s=0.5)
	    axes = plt.gca()
	    axes.set_xlim([0,gridpoints_x])
	    axes.set_ylim([0,gridpoints_y])
	plt.savefig(''.join([filename,'_pp_particle','.png']))

	return


def meshPlot( elements_x, elements_y, x_start, x_end, y_start, y_end, x_cluster, y_cluster, order ):

	# Plots the mesh of the simulation.

	# Defines size of grid.
        x = mp.mesh_generation(x_cluster,elements_x,x_start,x_end,order,2,'in')
        y = mp.mesh_generation(y_cluster,elements_y,y_start,y_end,order,2,'in')

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


def time_finder( filename, jump, total_timesteps ):

        total_files = total_timesteps/jump;
        range_vals = np.array(range(1,total_files))*jump
	time_old = 0

	for k in range_vals:
	
	    # Outputs counter to terminal.
	    print "Iteration no.: %s" % k

            # Reads data files.
            data,time,istep = rn.readnek(''.join([filename,'0.f',repr(k+1).zfill(5)]))

	    # Compute timestep.
	    dt = time - time_old
	    time_old = time

            print "Time step:      %s" % dt
	    print "Actual time:    %s" % time

	return

def linePlot( filename, order, start_file, jump, total_timesteps, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, particles ):
        # Plots data from a Nek5000 run.  Inputs are as follows:
        # filename: name that comes before the 0.f##### in the output files from Nek5000.
        # Order of the polynomial used for the spectral element method.
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
            data,time,istep = rn.readnek(''.join([filename,'0.f',repr(k).zfill(5)]))
            # Reshapes data onto uniform grid.
            [ mesh ] = rn.reshapenek(data, elements_y, elements_x)
            # Consider only the necessary number of plots.
            if (particles == 1):
#                temperature = mesh[:,:,3]
		horVel = mesh[:,:,0]
		verVel = mesh[:,:,1]
#		magVel = np.sqrt(np.square(verVel) + np.square(horVel))
		pressure = mesh[:,:,2]
            else:
#                temperature = mesh[:,:,5]
		horVel = mesh[:,:,2]
		verVel = mesh[:,:,3]
#		magVel = np.sqrt(np.square(verVel) + np.square(horVel))
		pressure = mesh[:,:,4]

            # Defines size of grid.
            x = mp.mesh_generation(x_cluster,elements_x,gridpoints_x,order,4,'out')
            y = mp.mesh_generation(y_cluster,elements_y,gridpoints_y,order,2,'out')

#	    x_data = temperature[200,:]
#	    name = '_tempLine_'
#	    x1 = 0
#	    x2 = 10
#	    y1 = 0
#	    y2 = 1
#	    pt.myPlot(x,x_data,time,'Horizontal Position','Temperature',filename,name,k/jump,x1,x2,y1,y2)

	    y_data = horVel[:,0]
	    name = '_wallHorVel_left_'
            x1 = -0.1
            x2 = 0.1
            y1 = 0
            y2 = 10
	    orientation = 'thin'
	    pt.myPlot(y_data,y,time,'Horizontal Velocity','Height',filename,name,k/jump,x1,x2,y1,y2,orientation)
            
	    x_data = verVel[-1,:]
            name = '_wallVerVel_top_'
            x1 = 0
            x2 = 10
            y1 = -1
            y2 = 1
	    orientation = 'long'
            pt.myPlot(x,x_data,time,'Horizontal Position','Vertical Velocity',filename,name,k/jump,x1,x2,y1,y2,orientation)

            x_data = pressure[-1,:]
            name = '_wallPress_top_'
            x1 = 0
            x2 = 10
            y1 = -0.05
            y2 = 0.05
	    orientation = 'long'
            pt.myPlot(x,x_data,time,'Horizontal Position','Pressure',filename,name,k/jump,x1,x2,y1,y2,orientation)

        return

