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

def PseudoColourPlotting( filename, order, start_file, jump, total_timesteps, numPlots, elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, gridType, particles ):
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
	    	if (numPlots == 1):
		    temperature = mesh[:,:,3]
	    	elif (numPlots == 2):
		    temperature = mesh[:,:,3]
		    verVel = mesh[:,:,1]
	    	else:
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
                else:
                    verVel = mesh[:,:,3]
                    horVel = mesh[:,:,2]
                    temperature = mesh[:,:,5]
                    magVel = np.sqrt(np.square(verVel) + np.square(horVel))

	    # Defines size of grid.
	    x = mp.mesh_generation(x_cluster,elements_x,gridpoints_x,order,4,'out')
	    y = mp.mesh_generation(y_cluster,elements_y,gridpoints_y,order,2,'out')
	    
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
		    dataPlot = [ abs(t)**0.5 for t in temperature ]
		    c_min = 0
		    c_max = 1
		    name = '_temp_'
		elif (plotNum == 1):
		    dataPlot = verVel
                    c_min = 0
                    c_max = 1
		    name = '_verVel_'
	    	elif (plotNum == 2):
                    dataPlot = magVel
                    c_min = 0
                    c_max = 1.0
		    name = '_magVel_'
                elif (plotNum == 3):
                    dataPlot = horVel
                    c_min = -0.2
                    c_max = 0.2
                    name = '_horVel_'

		if (particles == 0):
		    # Plots reshaped data
		    pt.myPcolour(np.transpose(x),y,dataPlot,time,'Horizontal position','Vertical position',range(0,11,10),range(0,11,10),filename,name,k/jump,vmin=c_min,vmax=c_max,cmap='RdBu_r')
		elif (particles == 1):
		    pt.particlePcolour(np.transpose(x),y,dataPlot,time,'Horizontal position','Vertical position',range(0,101,10),range(0,101,10),filename,name,k/jump,x_pos,y_pos,vmin=c_min,vmax=c_max)

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
		dx = points_x[x] - points_x[x-1]
		dy = points_y[y] - points_y[y-1]
		trap = sum1*dx*dy/4
		sum_total = sum_total + trap

	return sum_total


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
	buoyancy_flux = np.zeros(total_files+1)
	time_vec = np.zeros(total_files+1)
	all_energies = np.zeros(total_files+1)

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
	    magVel = np.sqrt(np.square(verVel) + np.square(horVel))
            temperature = mesh[:,:,3]
	    temperature = temperature + 273.15
            # Defines size of grid.
	    y_length,x_length = temperature.shape
	    
	    time_vec[file_counter] = time

	    # Defines size of grid.
            if( gridType == 0 ):
                [ x ] = geometricRatio(1/x_cluster,elements_x,gridpoints_x)
            else:
                [ x1 ] = geometricRatio(1/x_cluster,elements_x/2,gridpoints_x/2)
                [ x2 ] = geometricRatio(x_cluster,elements_x/2,gridpoints_x/2)
                x = np.concatenate([ x1[:-1], [ x+50 for x in x2 ] ])

	    [ y ] = geometricRatio(y_cluster,elements_y,gridpoints_y)

	    # Computing the integral of the energy.

	    density = [1027*(2+273)/T for T in temperature]
	    energy = np.square(magVel)
	    energy = 0.5*np.multiply(density,energy)
	    energy_at_time = trapezium(x,y,energy)

	    all_energies[file_counter] =  energy_at_time

	    file_counter = file_counter + 1

        vec = np.stack((time_vec, all_energies),axis=1)
        f = open("kinetic_energy.txt","w")
	f.write("# Time, Kinetic Energy`\n")
        np.savetxt(f, vec)

	return

 
def meshInd():

	for lamturb in range(0,2):

	    if (lamturb == 0):
    	        time_2040, ke_2040 = np.loadtxt("./save_lam_v1_20x40/kinetic_energy.txt", unpack = True)
                time_3060, ke_3060 = np.loadtxt("./save_lam_v2_30x60/kinetic_energy.txt", unpack = True)
	        time_4080, ke_4080 = np.loadtxt("./save_lam_v3_40x80/kinetic_energy.txt", unpack = True)
	    elif (lamturb == 1):
		time_2040, ke_2040 = np.loadtxt("./save_turb_v1_20x40/kinetic_energy.txt", unpack = True)
                time_3060, ke_3060 = np.loadtxt("./save_turb_v2_30x60/kinetic_energy.txt", unpack = True)
                time_4080, ke_4080 = np.loadtxt("./save_turb_v3_40x80/kinetic_energy.txt", unpack = True)

	    # Finds time such that all simulations have been run to the same (ish) time.
	    min_time =  min(max(time_2040),max(time_3060),max(time_4080))
	    length_2040 = tools.find_nearest(time_2040,min_time)
	    length_3060 = tools.find_nearest(time_3060,min_time)
	    length_4080 = tools.find_nearest(time_4080,min_time)

            plt.figure(figsize=(50, 30)) # Increases resolution.
	    ax = plt.axes()
	    plt.xlabel('Time',fontsize=80)
	    plt.ylabel('Kinetic Energy',fontsize=80)
	    plt.tick_params(axis='both', which='major', labelsize=60)
	    plt.tick_params(axis='both', which='minor', labelsize=60)
	    plt.plot(time_2040[2:length_2040],ke_2040[2:length_2040], label="Grid 20x40", linewidth = 5.0)
	    plt.plot(time_3060[2:length_3060],ke_3060[2:length_3060], label="Grid 30x60", linewidth = 5.0)
	    plt.plot(time_4080[2:length_4080],ke_4080[2:length_4080], label="Grid 40x80", linewidth = 5.0)
	    ax.yaxis.get_offset_text().set_fontsize(50)
	    plt.legend(fontsize=40)
	    if (lamturb == 0):
	        plt.savefig(''.join(['plume_v9_meshInd_keTimeCompare_laminar.png']))
	    elif (lamturb == 1):
		plt.savefig(''.join(['plume_v9_meshInd_keTimeCompare_turbulent.png']))
	    plt.close('all')

	    # Total kinetic energy of the system computed using the trapezium rule.
	    tot_ke_2040 = 0
	    tot_ke_3060 = 0
	    tot_ke_4080 = 0
	    for i in range(1,length_2040-1):
	        tot_ke_2040 = tot_ke_2040 + (time_2040[i+2]-time_2040[i+1])*(ke_2040[i+1]+ke_2040[i+2])/2
	    for i in range(1,length_3060-1):
                tot_ke_3060 = tot_ke_3060 + (time_3060[i+2]-time_3060[i+1])*(ke_3060[i+1]+ke_3060[i+2])/2
	    for i in range(1,length_4080-1):
                tot_ke_4080 = tot_ke_4080 + (time_4080[i+2]-time_4080[i+1])*(ke_4080[i+1]+ke_4080[i+2])/2	
	    elements = [800,1800,3200]
	    total_ke = [tot_ke_2040,tot_ke_3060,tot_ke_4080]

	    if (lamturb == 0):
                lam_ke = total_ke
            elif (lamturb == 1):
                turb_ke = total_ke

	    plt.figure(figsize=(50, 30)) # Increases resolution.
            ax = plt.axes()
            plt.xlabel('Number of Elements',fontsize=80)
            plt.ylabel('Total Kinetic Energy',fontsize=80)
            plt.tick_params(axis='both', which='major', labelsize=60)
            plt.tick_params(axis='both', which='minor', labelsize=60)
            plt.plot(elements,total_ke,linewidth = 5.0)
            ax.yaxis.get_offset_text().set_fontsize(50)
	    if (lamturb == 0):
                plt.savefig(''.join(['plume_v9_meshInd_keElementCompare_laminar.png']))
	    elif (lamturb == 1):
		plt.savefig(''.join(['plume_v9_meshInd_keElementCompare_turbulent.png']))
            plt.close('all')

	plt.figure(figsize=(50, 30)) # Increases resolution.
        ax = plt.axes()
        plt.xlabel('Number of Elements',fontsize=80)
        plt.ylabel('Total Kinetic Energy',fontsize=80)
        plt.tick_params(axis='both', which='major', labelsize=60)
        plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(elements,lam_ke,label="Laminar",linewidth = 5.0)
	plt.plot(elements,turb_ke,label="Turbulent",linewidth = 5.0)
	plt.legend(fontsize=40)
	#ax.set_xscale('log')
        ax.yaxis.get_offset_text().set_fontsize(50)
        plt.savefig(''.join(['plume_v9_meshInd_keElementCompare_both.png']))


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


def meshPlot( elements_x, elements_y, gridpoints_x, gridpoints_y, x_cluster, y_cluster, order ):

	# Plots the mesh of the simulation.

	# Defines size of grid.
        x = mp.mesh_generation(x_cluster,elements_x,gridpoints_x,order,4,'out')
        y = mp.mesh_generation(y_cluster,elements_y,gridpoints_y,order,2,'out')

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
