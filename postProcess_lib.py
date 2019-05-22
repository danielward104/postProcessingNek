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
import pickle

# Import user-defined modules.
import readingNek as rn
import mappings as mp
import plottingTools as pt
import generalTools as tools

def PseudoColourPlotting( filename, order, dimension, start_file, jump, final_timestep, numPlots, elements_x, elements_y, elements_z, y_slice, particles, simulation_currently_running):
        # Plots data from a Nek5000 run.  Inputs are as follows:
        # filename: name that comes before the 0.f##### in the output files from Nek5000.
        # Order of the polynomial used for the spectral element method.
        # start_file: file number to start at (usually leave at 1).
        # jump: number of 0.f##### files to skip between each plot.
        # final_timestep: number of the last 0.f##### file to consider.
        # numPlots: number of plots to produce (1 - temperature only, 2 - temperature and vertical velocity, 3 - temperature, vertical velocity, and magnitude of velocity).
        # elements_i: number of elements in the i-direction.
        # z_slice: location of the slice through the x-direction.
        # particles: switch to plot particle paths if particles are included in the simulation.

        order = order - 1

        applyMaxAndMinToPlots = 1

        quiver = 0      

        final_file = int(final_timestep/jump)

        #file_counter = 1
        
        if (start_file == 1):
            range_vals = [x - (jump - 1) for x in np.array(range(1,final_file+1))*jump]

        else:
            range_vals = np.array(range(int(math.floor(float(start_file)/jump)),final_file+1))*jump
            print("Make sure calculation of file numbers is correct")
            print(range_vals)

        # Reading in mesh data.
        if (simulation_currently_running == 0):
            data,time,istep,header,elmap,u_i,v_i,w_i,t_i,s_i = rn.readnek(''.join(['./dataFiles1/', \
                filename,'0.f00001']))
        else:
            data,time,istep,header,elmap,u_i,v_i,w_i,t_i,s_i = rn.readnek(''.join([filename,'0.f00001']))

        [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

        y = mesh[0,:,0,1]
        y_start = y[0]
        y_end = y[-1]

        # Find the point in the mesh where y_slice lies.
        nodes_y = elements_y*order + 1
        y_mesh = int(nodes_y*(1 - (y_end - y_slice)/(y_end - y_start)))

        # Defining x,y coordinates.
        x = mesh[:,y_mesh,0,0]
        z = mesh[0,y_mesh,:,2]

        x_start = x[0]
        x_end = x[-1]
        z_start = z[0]
        z_end = z[-1]

        # Plotting mesh.
        pt.meshPlot(x,z)

        for k in range_vals:

            file_num = int((k-1)/jump + 1)

            # Outputs counter to terminal.
            if (start_file == 1):
                files_remaining = int(final_file - k/jump)
            else:
                files_remaining = int(final_file - (k-range_vals[0])/jump - start_file/jump)

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            if (simulation_currently_running == 0):
                data,time,istep,header,elmap,u_i,v_i,w_i,t_i,s_i = rn.readnek(''.join(['./dataFiles1/', \
                    filename,'0.f',repr(k).zfill(5)]))
            else:
                data,time,istep,header,elmap,u_i,v_i,w_i,t_i,s_i = rn.readnek(''.join([filename, \
                    '0.f',repr(k).zfill(5)]))

            # Reshapes data onto uniform grid.
            if (dimension == 2):
                [ mesh ] = rn.reshapenek2D(data, elements_z, elements_x)
            
            elif (dimension == 3):

                [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)
                mesh = mesh[:,y_mesh,:,:]
                
            # For some reason, coordinates are not plotted in files following the first when particles 
            # are included in the simulation.
            if (particles == 1):
                if (k > 1):
                    u_i = u_i - 3
                    v_i = v_i - 3
                    w_i = w_i - 3
                    t_i = t_i - 3

            # Consider only the necessary number of plots.
            if (numPlots == 1):
                temperature = np.transpose(mesh[:,:,t_i])
                scalar = np.transpose(mesh[:,:,s_i])

            elif (numPlots == 2):
                temperature = np.transpose(mesh[:,:,t_i])
                scalar = np.transpose(mesh[:,:,s_i])
                horVel = np.transpose(mesh[:,:,u_i])
                horVel2 = np.transpose(mesh[:,:,v_i])
                verVel = np.transpose(mesh[:,:,w_i])
                magVel = np.sqrt(np.square(verVel) + np.square(horVel) + np.square(horVel2))

######################################################################
############################ Particles!! #############################
######################################################################

            # Reading in particle data.
            if (particles == 1):
                npart = (k)*4
                pname = ''.join(['part',repr(npart).zfill(5),'.3D'])
               
                text_file = open(pname,'rb')                
                lines = text_file.read().decode()
                text_file.close()

                lines = lines.splitlines()

                x_pos = np.zeros(len(lines))
                y_pos = np.zeros(len(lines))
                z_pos = np.zeros(len(lines))

                for i in range(0,len(lines)):
                    line = lines[i].split()
                    x_pos[i] = float(line[0])
                    y_pos[i] = float(line[1])
                    z_pos[i] = float(line[2])

                dataPlot = temperature 
#                c_min = 0.0
#                c_max = 1.0
                name = 'temperature'
                pt.particlePcolour(np.transpose(x),z,dataPlot,time,'Horizontal position', \
                        'Vertical position',filename,name,file_num,x_pos,z_pos,cmap='RdBu_r') #, \
#                        vmin=c_min,vmax=c_max)
            
#                pt.myPcolour(np.transpose(x),y,np.transpose(dataPlot),time,'Horizontal position', \
#                        'Vertical position',filename,name,file_num,cmap='RdBu_r', \
#                        vmin=c_min,vmax=c_max)
#
#
                dataPlot = verVel
#                c_min = 0.0
#                c_max = 0.1
                name = 'verVel'
                pt.particlePcolour(np.transpose(x),z,dataPlot,time,'Horizontal position', \
                        'Vertical position',filename,name,file_num,x_pos,z_pos,cmap='RdBu_r') #, \
#                        vmin=c_min,vmax=c_max)
        
                dataPlot = horVel
                name = 'horVel'
                pt.particlePcolour(np.transpose(x),z,dataPlot,time,'Horizontal position', \
                        'Vertical position',filename,name,file_num,x_pos,z_pos,cmap='RdBu_r')

                dataPlot = magVel
                name = 'magVel'
                pt.particlePcolour(np.transpose(x),z,dataPlot,time,'Horizontal position', \
                        'Vertical position',filename,name,file_num,x_pos,z_pos,cmap='RdBu_r') 

######################################################################
########################## No particles!! ############################
######################################################################

            else:
                for plotNum in range(0,numPlots):

                
                    if (plotNum == 0):

                        dataPlot = scalar
                        name = 'passive-scalar'

                        x_plot = x
                        z_plot = z

                        x_plot_start = x_start
                        x_plot_end = x_end
                        z_plot_start = z_start
                        z_plot_end = z_end

                        zoomify = 0

                        if (zoomify == 1):
                                
                            x_i_start = 0
                            x_i_end = len(x)
                            z_i_start = 0
                            z_i_end = int(len(z)/1.5)

                            x_plot = x[x_i_start:x_i_end]
                            z_plot = z[z_i_start:z_i_end]

                            x_plot_start = x_plot[0]
                            x_plot_end = x_plot[-1]
                            z_plot_start = z_plot[0]
                            z_plot_end = z_plot[-1]
                                
                            dataPlot = np.array(dataPlot)
                            dataPlot = dataPlot[z_i_start:z_i_end,x_i_start:x_i_end]
                        
                            name = 'temperature_zoom'

                        c_min = 0
                        c_max = 1

                        if (applyMaxAndMinToPlots == 1):
                                pt.myPcolour(np.transpose(x_plot),z_plot,dataPlot,time,\
                                x_plot_start,x_plot_end,z_plot_start,z_plot_end,\
                                'Horizontal position','Vertical position',filename,name,\
                                file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)

                        else:
                                pt.myPcolour(np.transpose(x_plot),z_plot,dataPlot,time,\
                                x_plot_start,x_plot_end,z_plot_start,z_plot_end,\
                                'Horizontal position','Vertical position',filename,name,\
                                file_num,cmap='RdBu_r')

#                        pt.particlePcolour(np.transpose(x),y,dataPlot,time,'Horizontal position', \
#                                'Vertical position',filename,name,file_num,x_pos,y_pos, \
#                                vmin=c_min,vmax=c_max,cmap='RdBu_r')

#                       if (quiver == 1):
#                       quiverPlotx = horVel
#                       quiverPloty = verVel 

                    elif (plotNum == 1):

                        # Plotting Vertical velocity
                        dataPlot = verVel
                        c_min = -2.0
                        c_max = 8.0

                        #bound = np.amax(abs(verVel))
                        #c_min = -bound
                        #c_max = bound

                        name = 'z-velocity'

                        x_plot = x
                        z_plot = z

                        x_plot_start = x_start
                        x_plot_end = x_end
                        z_plot_start = z_start
                        z_plot_end = z_end

                        if (zoomify == 1):

                            x_i_start = 0
                            x_i_end = len(x)
                            z_i_start = 0
                            z_i_end = int(len(z)/1.5)

                            x_plot = x[x_i_start:x_i_end]
                            z_plot = z[z_i_start:z_i_end]

                            x_plot_start = x_plot[0]
                            x_plot_end = x_plot[-1]
                            z_plot_start = z_plot[0]
                            z_plot_end = z_plot[-1]

                            dataPlot = np.array(dataPlot)
                            dataPlot = dataPlot[z_i_start:z_i_end,x_i_start:x_i_end]

                            name = 'z-velocity_zoom'


                        if (applyMaxAndMinToPlots == 1):
                                pt.myPcolour(np.transpose(x_plot),z_plot,dataPlot,time,\
                                x_plot_start,x_plot_end,z_plot_start,z_plot_end,\
                                'Horizontal position','Vertical position',filename,name,\
                                file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)

                        else:
                                pt.myPcolour(np.transpose(x_plot),z_plot,dataPlot,time,\
                                x_plot_start,x_plot_end,z_plot_start,z_plot_end,\
                                'Horizontal position','Vertical position',filename,name,\
                                file_num,cmap='RdBu_r')

                        
                        # Plotting magnitude of velocity
                        #dataPlot = magVel
                        #c_min = 0
                        #c_max = 1
                        #name = 'velocity-magnitude'

                        #if (zoomify == 1):
                        #    name = 'velocity-magnitude_zoom'
#
#                            dataPlot = np.array(dataPlot)
#                            dataPlot = dataPlot[y_i_start:y_i_end,x_i_start:x_i_end]
#
#                        pt.myPcolour(np.transpose(x_plot),y_plot,dataPlot,time,\
#                                x_plot_start,x_plot_end,y_plot_start,y_plot_end \
#                                ,'Horizontal position','Vertical position',filename,name \
#                                ,file_num,cmap='RdBu_r')#,vmin=c_min,vmax=c_max)
#
                        # Plotting horizontal velocity
                        dataPlot = horVel
                        c_min = -0.1
                        c_max = 0.1
                        name = 'x-velocity'

                        if (zoomify == 1):
                            name = 'x-velocity_zoom'

                            dataPlot = np.array(dataPlot)
                            dataPlot = dataPlot[z_i_start:z_i_end,x_i_start:x_i_end]

                        pt.myPcolour(np.transpose(x_plot),z_plot,dataPlot,time,\
                                x_plot_start,x_plot_end,z_plot_start,z_plot_end \
                                ,'Horizontal position','Vertical position',filename,name \
                                ,file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)

        return


def average_field( filename, order, dimension, start_file, jump, final_timestep, elements_x, elements_y, elements_z ):

        final_file = int(final_timestep/jump)

        if (start_file == 1):
            range_vals = [x - (jump - 1) for x in np.array(range(1,final_file+1))*jump]
        else:
            range_vals = np.array(range(int(math.floor(float(start_file)/jump)),final_file+1))*jump
            print("Make sure calculation of file numbers is correct")
            print(range_vals)

        # Reading in mesh data.
        data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join([filename,'0.f00001']))
        [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

        # Defining x,y,z coordinates.
        x = mesh[:,0,0,0]
        y = mesh[0,:,0,1]
        z = mesh[0,0,:,2]
        
        x_start = x[0]
        x_end = x[-1]
        y_start = y[0]
        y_end = y[-1]
        z_start = z[0]
        z_end = z[-1]

        xlength = x.shape
        xlength = xlength[0]
        x00 = int((xlength - 1)/2)
        zlength = z.shape
        zlength = zlength[0]
        z00 = int((zlength - 1)/2)
        ylength = y.shape
        ylength = ylength[0]

        u_r_avg = np.zeros((ylength,x00+1))
        u_y_avg = np.zeros((ylength,x00+1))
        counter = 0

        for k in range_vals:

            file_num = int((k-1)/jump + 1)

            # Outputs counter to terminal.
            if (start_file == 1):
                files_remaining = int(final_file - k/jump)
            else:
                files_remaining = int(final_file - (k-range_vals[0])/jump - start_file/jump)

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join([filename,'0.f',repr(k).zfill(5)]))

            # Reshapes data onto uniform grid.
            if (dimension == 2):
                [ mesh ] = rn.reshapenek2D(data, elements_y, elements_x)

            elif (dimension == 3):
                [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

            u_x = np.transpose(mesh[:,:,:,u_i])
            u_y = np.transpose(mesh[:,:,:,v_i])
            u_z = np.transpose(mesh[:,:,:,w_i])

            #u_r = np.divide(np.multiply(x,u_x) + np.multiply(y,u_y),np.sqrt(np.multiply(x,x) + np.multiply(y,y)))

            xlength = x.shape
            xlength = xlength[0]
            x00 = int((xlength - 1)/2)
            zlength = z.shape
            zlength = zlength[0]
            z00 = int((zlength - 1)/2)
            
            u_r_1 =                         u_x[x00,        :,z00:zlength]
            u_r_2 = np.fliplr(             -u_x[x00,        :,0:z00+1])
            u_r_3 = np.transpose(           u_z[x00:xlength,:,z00])
            u_r_4 = np.transpose(np.flipud(-u_z[0:x00+1,    :,z00]))
 
            u_y_1 =                        u_y[x00,        :,z00:zlength]
            u_y_2 = np.fliplr(             u_y[x00,        :,0:z00+1])
            u_y_3 = np.transpose(          u_y[x00:xlength,:,z00])
            u_y_4 = np.transpose(np.flipud(u_y[0:x00+1,    :,z00]))

            u_r_avg = u_r_avg + (u_r_1 + u_r_2 + u_r_3 + u_r_4)/4
            u_y_avg = u_y_avg + (u_y_1 + u_y_2 + u_y_3 + u_y_4)/4
            counter = counter + 1

#            c_min = -0.01
#            c_max = 0.01
#
#            name = 'u_r_1'
#            pt.myPcolour(np.transpose(x[x00:xlength]),y,u_r_1,time,\
#                    x[x00],x[-1],y[0],y[-1],\
#                    'Horizontal position','Vertical position',filename,name,\
#                    file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)
#
#            name = 'u_r_2'
#            pt.myPcolour(np.transpose(x[x00:xlength]),y,u_r_2,time,\
#                    x[x00],x[-1],y[0],y[-1],\
#                    'Horizontal position','Vertical position',filename,name,\
#                    file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)
#
#            name = 'u_r_3'
#            pt.myPcolour(np.transpose(x[x00:xlength]),y,u_r_3,time,\
#                    x[x00],x[-1],y[0],y[-1],\
#                    'Horizontal position','Vertical position',filename,name,\
#                    file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)
#
#            name = 'u_r_4'
#            pt.myPcolour(np.transpose(x[x00:xlength]),y,u_r_4,time,\
#                    x[x00],x[-1],y[0],y[-1],\
#                    'Horizontal position','Vertical position',filename,name,\
#                    file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)
#
#            name = 'u_r_avg'
#            pt.myPcolour(np.transpose(x[x00:xlength]),y,u_r_avg,time,\
#                    x[x00],x[-1],y[0],y[-1],\
#                    'Horizontal position','Vertical position',filename,name,\
#                    file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)


        u_r_avg = u_r_avg/counter
        u_y_avg = u_y_avg/counter
        c_min = -0.01
        c_max = 0.01

        name = 'u_r_avg'
        pt.myPcolour(np.transpose(x[x00:xlength]),y,u_r_avg,time,\
                x[x00],x[-1],y[0],y[-1],\
                'Horizontal position','Vertical position',filename,name,\
                file_num,cmap='RdBu_r',vmin=c_min,vmax=c_max)

        
        name = 'u_y_avg'
        pt.myPcolour(np.transpose(x[x00:xlength]),y,u_y_avg,time,\
                x[x00],x[-1],y[0],y[-1],\
                'Horizontal position','Vertical position',filename,name,\
                file_num,cmap='RdBu_r')#,vmin=c_min,vmax=c_max)

        f = open('avg_r_vel.file','wb')
        pickle.dump(u_r_avg,f)
        f.close()

        f = open('avg_y_vel.file','wb')
        pickle.dump(u_r_avg,f)
        f.close()


        return

def TKE( filename, order, dimension, start_file, jump, final_timestep, elements_x, elements_y, elements_z, simulation_currently_running):

        final_file = int(final_timestep/jump)

        u_r_avg = pickle.load(open('avg_r_vel.file','rb'))
        u_y_avg = pickle.load(open('avg_y_vel.file','rb'))

        if (start_file == 1):
            range_vals = [x - (jump - 1) for x in np.array(range(1,final_file+1))*jump]
        else:
            range_vals = np.array(range(int(math.floor(float(start_file)/jump)),final_file+1))*jump
            print("Make sure calculation of file numbers is correct")
            print(range_vals)

        # Reading in mesh data.
        if (simulation_currently_running == 0):
            data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join(['./data_files/', \
                filename,'0.f00001']))
        else:
            data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join([filename,'0.f00001']))

        [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

        # Defining x,y,z coordinates.
        x = mesh[:,0,0,0]
        y = mesh[0,:,0,1]
        z = mesh[0,0,:,2]
        
        x_start = x[0]
        x_end = x[-1]
        y_start = y[0]
        y_end = y[-1]
        z_start = z[0]
        z_end = z[-1]

        xlength = x.shape
        xlength = xlength[0]
        x00 = int((xlength - 1)/2)
        zlength = z.shape
        zlength = zlength[0]
        z00 = int((zlength - 1)/2)
        ylength = y.shape
        ylength = ylength[0]

        u_r_avg = np.zeros((ylength,x00+1))
        counter = 0

        for k in range_vals:

            file_num = int((k-1)/jump + 1)

            # Outputs counter to terminal.
            if (start_file == 1):
                files_remaining = int(final_file - k/jump)
            else:
                files_remaining = int(final_file - (k-range_vals[0])/jump - start_file/jump)

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            if (simulation_currently_running == 0):
                data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join(['./data_files/', \
                    filename,'0.f',repr(k).zfill(5)]))
            else:
                data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join([filename, \
                    '0.f',repr(k).zfill(5)]))

            # Reshapes data onto uniform grid.
            if (dimension == 2):
                [ mesh ] = rn.reshapenek2D(data, elements_y, elements_x)

            elif (dimension == 3):
                [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)
            

            u_x = np.transpose(mesh[:,:,:,u_i])
            u_y = np.transpose(mesh[:,:,:,v_i])
            u_z = np.transpose(mesh[:,:,:,w_i])

            xlength = x.shape
            xlength = xlength[0]
            x00 = int((xlength - 1)/2)
            zlength = z.shape
            zlength = zlength[0]
            z00 = int((zlength - 1)/2)
            
            u_r_1 =                         u_x[x00,        :,z00:zlength]
            u_r_2 = np.fliplr(             -u_x[x00,        :,0:z00+1])
            u_r_3 = np.transpose(           u_z[x00:xlength,:,z00])
            u_r_4 = np.transpose(np.flipud(-u_z[0:x00+1,    :,z00]))
 
            u_y_1 =                        u_y[x00,        :,z00:zlength]
            u_y_2 = np.fliplr(             u_y[x00,        :,0:z00+1])
            u_y_3 = np.transpose(          u_y[x00:xlength,:,z00])
            u_y_4 = np.transpose(np.flipud(u_y[0:x00+1,    :,z00]))
            
            u_r_prime = abs(u_r_avg - u_r_1)
            u_y_prime = abs(u_y_avg - u_y_1)

            k = 0.5*(np.square(u_r_prime) + np.square(u_y_prime))
            k = trapezium_2D( x[x00:xlength], y, k )

            
        

        return


def trapezium_1D( points_x, data ):

        # Computes the approximate two-dimensional integral of the function represented by 'data'.

        sum_total = 0
        x_tot = np.shape(points_x)
        x_tot = x_tot[0] - 1

        for x in range(0,x_tot):
            sum1 = data[x] + data[x+1]
            dx = points_x[x+1] - points_x[x]
            trap = sum1*dx/2
            sum_total = sum_total + trap

        return sum_total



def trapezium_2D( points_x, points_y, data ):

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

def trapezium_3D( points_x, points_y, points_z, data ):

        # Computes the approximate three-dimensional integral of the function represented by 'data'.
    
        sum_total = 0
        x_tot = np.shape(points_x)
        x_tot = x_tot[0] - 1
        y_tot = np.shape(points_y)
        y_tot = y_tot[0] - 1
        z_tot = np.shape(points_z)
        z_tot = z_tot[0] - 1

        for x in range(0,x_tot):
            for y in range(0,y_tot):
                for z in range(0,z_tot):
                    sum1 = data[x,y,z] + data[x+1,y,z] + data[x,y+1,z] + data[x+1,y+1,z] + \
                            data[x,y,z+1] + data[x+1,y,z+1] + data[x,y+1,z+1] + data[x+1,y+1,z+1]

                    dx = points_x.item(x+1) - points_x.item(x)
                    dy = points_y.item(y+1) - points_y.item(y)
                    dz = points_z.item(z+1) - points_z.item(z)
                    trap = sum1*dx*dy*dz/8
                    sum_total = sum_total + trap
                    
        return sum_total


def integrateDomain( filename, order, dimension, jump, final_timestep, elements_x, elements_y, elements_z, x_start, x_end, y_start, y_end, z_start, z_end, x_cluster, y_cluster, gridType ):

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

        final_file = int(final_timestep/jump)
        ambientTemp = 273.15
        ambientDensity = 1027
        g = 9.81
        expCoeff = 0.0002

        range_vals = np.array(range(0,final_file))*jump

        # Initialises files to write to.  Column 1 will contain time data, column 2 will contain 
            # the data represented by the name of the file.
        f1 = open("kinetic_energy.txt","wb")
#        f2 = open("buoyancy.txt","wb")
#        f3 = open("avgVels.txt","wb")

        for k in range_vals:
            
            file_num = k/jump + 1
            
            # Outputs counter to terminal.
            files_remaining = int(final_file - k/jump)
            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            data,time,istep,header,elmap = rn.readnek(''.join([filename,'0.f',repr(k+1).zfill(5)]))
            # Reshapes data onto uniform grid.

            # Defines size of grid.
            x = mp.mesh_generation(x_cluster,elements_x,x_start,x_end,order,2,'in')
            y = mp.mesh_generation(y_cluster,elements_y,y_start,y_end,order,2,'in')

            if (dimension == 2):

                [ mesh ] = rn.reshapenek2D(data, elements_y, elements_x)

                verVel = mesh[:,:,3]
                horVel = mesh[:,:,2]
                magVel = np.sqrt(np.square(verVel) + np.square(horVel))
                temperature = mesh[:,:,5]
                temperature = temperature + 273.15

            elif (dimension == 3):

                [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)
                
                verVel = mesh[:,:,4,4]
                horVel = mesh[:,:,4,3]
                magVel = np.sqrt(np.square(verVel) + np.square(horVel))
                temperature = mesh[:,:,4,7]
                temperature = temperature + 273.15
                
                z = np.linspace(z_start,z_end,order*elements_z+1)

            # Computing the integral of the energy and buoyancy.

            density = [1027*(2+273)/T for T in temperature]
            density = [ambientDensity*(1-expCoeff*(T - ambientTemp)) for T in temperature]
            energy = np.square(magVel)
            energy = 0.5*np.multiply(density,energy)
            
            buoyancy = [g*(p-ambientDensity)/p for p in density]

            if (dimension == 3):
                energy_at_time = trapezium_2D(y,x,energy)
                buoyancy_total = trapezium_2D(y,x,np.array(buoyancy))
                avgVel = trapezium_2D(y,x,magVel)

            elif (dimension == 2):
                energy_at_time = trapezium_3D(x,y,z,energy)
                buoyancy_total = trapezium_3D(x,y,z,np.array(buoyancy))
                avgVel = trapezium_3D(x,y,z,magVel)

            # Writing data to file.
            np.savetxt(f1, np.atleast_2d(np.array([time,energy_at_time])), fmt='%1.8e')
#            np.savetxt(f2, np.atleast_2d(np.array([time,buoyancy_total])), fmt='%1.8e')
#            np.savetxt(f3, np.atleast_2d(np.array([time,avgVel])), fmt='%1.8e')

        # Closing files.
        f1.close()
        f2.close()
        f3.close()

        return


def integratePlume( filename, file_loc, order, dimension, start_file, jump, final_timestep, numPlots, elements_x, elements_y, elements_z ):

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

        order = order - 1

        final_file = int(final_timestep/jump)

        ii = 0

        ambientTemp = 273.15
        ambientDensity = 1027
        g = 9.81
        expCoeff = 0.0002

        if (start_file == 1):
            range_vals = [x - (jump - 1) for x in np.array(range(1,final_file+1))*jump]

        else:
            range_vals = np.array(range(int(math.floor(float(start_file)/jump)),final_file+1))*jump
            print("Make sure calculation of file numbers is correct")
            print(range_vals)

        # Reading in mesh data.
        data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join([file_loc, \
                filename,'0.f00001']))

        [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

        # Defining x,y,z coordinates.
        x = mesh[:,0,0,0]
        y = mesh[0,:,0,1]
        z = mesh[0,0,:,2]

        x_start = x[0]
        x_end = x[-1]
        y_start = y[0]
        y_end = y[-1]
        z_start = z[0]
        z_end = z[-1]

        height = len(z)

        x_len = len(x)
        y_len = len(y)

        radius = np.zeros((len(x),len(y)))

        for i in range(x_len):
            for j in range(y_len):
                radius[i,j] = (x[i]**2 + y[j]**2)**0.5

        # Initialises files to write to.  Column 1 will contain time data, column 2 will contain 
            # the data represented by the name of the file.
#        f1 = open("volume_flux.txt","wb")
#        f2 = open("buoyancy.txt","wb")
#        f3 = open("avgVels.txt","wb")

        for t in range_vals:

            file_num = int((t-1)/jump + 1)

            # Outputs counter to terminal.
            if (start_file == 1):
                files_remaining = int(final_file - t/jump)
            else:
                files_remaining = int(final_file - (t-range_vals[0])/jump - start_file/jump)

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            data,time,istep,header,elmap,u_i,v_i,w_i,t_i = \
                        rn.readnek(''.join([file_loc,filename,'0.f',repr(t).zfill(5)]))

            # Reshapes data onto uniform grid.
            if (dimension == 2):
                [ mesh ] = rn.reshapenek2D(data, elements_z, elements_x)

            elif (dimension == 3):

                [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

            temperature = mesh[:,:,:,t_i]
#            horVel = np.transpose(mesh[:,:,:,u_i])
            verVel = mesh[:,:,:,v_i]
#            magVel = np.sqrt(np.square(verVel) + np.square(horVel))

            Q = np.zeros((height,1))
            M = np.zeros((height,1))
            F = np.zeros((height,1))

            x_loc3D = []
            y_loc3D = []
            z_loc3D = []

            # Running through heights
            for k in range(3,height):

                int1 = np.multiply(verVel[:,:,k],radius)
                int2 = np.multiply(np.multiply(verVel[:,:,k],verVel[:,:,k]),radius)
                int3 = np.multiply(np.multiply(-temperature[:,:,k],verVel[:,:,k]),radius)
        
                sum_total1 = 0
                sum_total2 = 0
                sum_total3 = 0

                x_tot = x_len - 1
                y_tot = y_len - 1

                tol = 0.1

                for i in range(0,x_tot):
                    for j in range(0,y_tot):
                        if temperature[i,j,k] > tol:
                            sum1 = int1[i,j] + int1[i+1,j] + int1[i,j+1] + int1[i+1,j+1]
                            sum2 = int2[i,j] + int2[i+1,j] + int2[i,j+1] + int2[i+1,j+1]
                            sum3 = int3[i,j] + int3[i+1,j] + int3[i,j+1] + int3[i+1,j+1]
                            dx = x[i+1] - x[i]
                            dy = y[j+1] - y[j]
                            trap1 = sum1*dx*dy/4
                            trap2 = sum2*dx*dy/4
                            trap3 = sum3*dx*dy/4
                            sum_total1 = sum_total1 + trap1
                            sum_total2 = sum_total2 + trap2
                            sum_total3 = sum_total3 + trap3

                Q[k] = sum_total1
                M[k] = sum_total2
                F[k] = sum_total3

            plt.figure(figsize=(20,20))
            plt.plot(Q,z,color='blue')
            plt.plot(M,z,color='green')
            plt.plot(F,z,color='red')

            plt.gca().legend(('Q','M','F'))

            plt.savefig(''.join(['fluxes',repr(file_num).zfill(5),'.png']),bbox_inches='tight')
            plt.close('all')

            # Writing data to file.
            f = open(''.join(['./integration/volumeFlux/volumeFlux',repr(file_num).zfill(5),'.txt']),\
                                                                                                "wb")
            np.savetxt(f, Q, fmt='%1.8e', header=''.join(['time = ',repr(time)]))
            f.close()


            # Computing the integral of the energy and buoyancy.

#            density = [1027*(2+273)/T for T in temperature]
#            density = [ambientDensity*(1-expCoeff*(T - ambientTemp)) for T in temperature]
#            energy = np.square(magVel)
#            energy = 0.5*np.multiply(density,energy)
#
#            buoyancy = [g*(p-ambientDensity)/p for p in density]

#            energy_at_time = trapezium_2D(y,x,energy)
#            buoyancy_total = trapezium_2D(y,x,np.array(buoyancy))
#            avgVel = trapezium_2D(y,x,magVel)

            # Computing volume flux.
#           radius = 
#           wr = verVel*radius

        return



def plotparticlesonly( order, start_file, jump, final_timestep, elements_x, elements_y, elements_z, x_start, x_end, y_start, y_end, z_slice, x_cluster, y_cluster ):


        final_file = final_timestep/jump;

        if (start_file == 1):
            range_vals = [x - (jump - 1) for x in np.array(range(1,final_file+1))*jump]

        else:
            range_vals = [x + 1 for x in np.array(range(int(math.floor(float(start_file)/jump)),\
                final_file+1))*jump]
            print("Make sure calculation of file numbers is correct")
            print(range_vals)

        for k in range_vals:

            file_num = k/jump + 1

            # Outputs counter to terminal.
            if (start_file == 1):
                files_remaining = final_file - k/jump
            else:
                files_remaining = final_file - (k-range_vals[0])/jump - start_file/jump

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Defines size of grid.
            x = mp.mesh_generation(x_cluster,elements_x,x_start,x_end,order,2,'in')
            y = mp.mesh_generation(y_cluster,elements_y,y_start,y_end,order,2,'in')

            # Reading in particle data.
            npart = k
            pname = ''.join(['part',repr(npart).zfill(5),'.3D'])
            text_file = open(pname,"r")

            lines = text_file.read().strip()
            lines = lines.splitlines()
            lines = np.asarray([ map(float, line.split()) for line in lines ])
            x_pos = lines[:,0]
            y_pos = lines[:,1]
            z_pos = lines[:,2]

            pt.particleOnlyPlot(x,y,'Horizontal position','Vertical position',\
                                        file_num,x_pos,y_pos)



        return
 
def meshInd():

        time_2010, ke_2010 = np.loadtxt("./LES_20x10x05/kinetic_energy.txt", unpack = True)
        time_4020, ke_4020 = np.loadtxt("./LES_40x20x10/kinetic_energy.txt", unpack = True)
        time_8040, ke_8040 = np.loadtxt("./LES_80x40x20/kinetic_energy.txt", unpack = True)
#        time_12060, ke_12060 = np.loadtxt("./save_v3_120x60/kinetic_energy.txt", unpack = True)
#        time_16080, ke_16080 = np.loadtxt("./save_v4_160x80/kinetic_energy.txt", unpack = True)

        time_2010, buoy_2010 = np.loadtxt("./LES_20x10x05/buoyancy.txt", unpack = True)
        time_4020, buoy_4020 = np.loadtxt("./LES_40x20x10/buoyancy.txt", unpack = True)
        time_8040, buoy_8040 = np.loadtxt("./LES_80x40x20/buoyancy.txt", unpack = True)
#        time_12060, buoy_12060 = np.loadtxt("./save_v3_120x60/buoyancy.txt", unpack = True)
#        time_16080, buoy_16080 = np.loadtxt("./save_v4_160x80/buoyancy.txt", unpack = True)

        time_2010, avgVel_2010 = np.loadtxt("./LES_20x10x05/avgVels.txt", unpack = True)
        time_4020, avgVel_4020 = np.loadtxt("./LES_40x20x10/avgVels.txt", unpack = True)
        time_8040, avgVel_8040 = np.loadtxt("./LES_80x40x20/avgVels.txt", unpack = True)
#        time_12060, avgVel_12060 = np.loadtxt("./save_v3_120x60/avgVels.txt", unpack = True)
#        time_16080, avgVel_16080 = np.loadtxt("./save_v4_160x80/avgVels.txt", unpack = True)

        # Finds time such that all simulations have been run to the same (ish) time.
        min_time =  min(max(time_2010),max(time_4020),max(time_8040))#,max(time_12060),max(time_16080))
        length_2010 = tools.find_nearest(time_2010,min_time)
        length_4020 = tools.find_nearest(time_4020,min_time)
        length_8040 = tools.find_nearest(time_8040,min_time)
#        length_12060 = tools.find_nearest(time_12060,min_time)
#        length_16080 = tools.find_nearest(time_16080,min_time)

        # Plots kinetic energy vs. time.
        plt.figure(figsize=(50, 30)) # Increases resolution.
        ax = plt.axes()
        plt.xlabel('Time',fontsize=80)
        plt.ylabel('Kinetic Energy',fontsize=80)
        plt.tick_params(axis='both', which='major', labelsize=60)
        plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(time_2010[2:length_2010],ke_2010[2:length_2010], label="Grid 20x10", linewidth = 5.0)
        plt.plot(time_4020[2:length_4020],ke_4020[2:length_4020], label="Grid 40x20", linewidth = 5.0)
        plt.plot(time_8040[2:length_8040],ke_8040[2:length_8040], label="Grid 80x40", linewidth = 5.0)
#        plt.plot(time_12060[2:length_12060],ke_12060[2:length_12060], label="Grid 120x60", linewidth = 5.0)
#        plt.plot(time_16080[2:length_16080],ke_16080[2:length_16080], label="Grid 160x80", linewidth = 5.0)

 #       plt.plot(time_4020,ke_4020, label="Grid 40x20", linewidth = 5.0)

        ax.yaxis.get_offset_text().set_fontsize(50)
        plt.legend(fontsize=40)
        plt.savefig(''.join(['plume_v3_stratified_keTime.png']),bbox_inches='tight')    
        plt.close('all')

        # Plots buoyancy vs. time.
        plt.figure(figsize=(50, 30)) # Increases resolution.
        ax = plt.axes()
        plt.xlabel('Time',fontsize=80)
        plt.ylabel('Buoyancy',fontsize=80)
        plt.tick_params(axis='both', which='major', labelsize=60)
        plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(time_2010[2:length_2010],buoy_2010[2:length_2010], label="Grid 20x10", linewidth = 5.0)
        plt.plot(time_4020[2:length_4020],buoy_4020[2:length_4020], label="Grid 40x20", linewidth = 5.0)
        plt.plot(time_8040[2:length_8040],buoy_8040[2:length_8040], label="Grid 80x40", linewidth = 5.0)
#        plt.plot(time_12060[2:length_12060],buoy_12060[2:length_12060], label="Grid 120x60", linewidth = 5.0)
#        plt.plot(time_16080[2:length_16080],buoy_16080[2:length_16080], label="Grid 160x80", linewidth = 5.0)
        ax.yaxis.get_offset_text().set_fontsize(50)
        plt.legend(fontsize=40)
        plt.savefig(''.join(['plume_v3_stratified_buoyTime.png']),bbox_inches='tight')
        plt.close('all')

        # Plots average velocity.

        plt.figure(figsize=(50, 30)) # Increases resolution.
        ax = plt.axes()
        plt.xlabel('Time',fontsize=80)
        plt.ylabel('Average Velocity',fontsize=80)
        plt.tick_params(axis='both', which='major', labelsize=60)
        plt.tick_params(axis='both', which='minor', labelsize=60)
        plt.plot(time_2010[2:length_2010],avgVel_2010[2:length_2010], label="Grid 20x10", linewidth = 5.0)
        plt.plot(time_4020[2:length_4020],avgVel_4020[2:length_4020], label="Grid 40x20", linewidth = 5.0)
        plt.plot(time_8040[2:length_8040],avgVel_8040[2:length_8040], label="Grid 80x40", linewidth = 5.0)
#        plt.plot(time_12060[2:length_12060],avgVel_12060[2:length_12060], label="Grid 120x60", linewidth = 5.0)
#        plt.plot(time_16080[2:length_16080],avgVel_16080[2:length_16080], label="Grid 160x80", linewidth = 5.0)
        ax.yaxis.get_offset_text().set_fontsize(50)
        plt.legend(fontsize=40)
        plt.savefig(''.join(['plume_v3_stratified_avgVel.png']),bbox_inches='tight')
        plt.close('all')


        # Total kinetic energy of the system computed using the trapezium rule.
#        tot_ke_4020 = 0
#        tot_ke_8040 = 0
#        tot_ke_12060 = 0
#        for i in range(1,length_4020-1):
#            tot_ke_4020 = tot_ke_4020 + (time_4020[i+2]-time_4020[i+1])*(ke_4020[i+1]+ke_4020[i+2])/2
#        for i in range(1,length_8040-1):
#            tot_ke_8040 = tot_ke_8040 + (time_8040[i+2]-time_8040[i+1])*(ke_8040[i+1]+ke_8040[i+2])/2
#        for i in range(1,length_12060-1):
#            tot_ke_12060 = tot_ke_12060 + (time_12060[i+2]-time_12060[i+1])*(ke_12060[i+1]+ke_12060[i+2])/2     

#        total_ke = [tot_ke_4020,tot_ke_8040,tot_ke_12060]

#        plt.figure(figsize=(50, 30)) # Increases resolution.
#        ax = plt.axes()
#        plt.xlabel('Number of Elements',fontsize=80)
#        plt.ylabel('Total Kinetic Energy',fontsize=80)
#        plt.tick_params(axis='both', which='major', labelsize=60)
#        plt.tick_params(axis='both', which='minor', labelsize=60)
#        plt.plot(elements,total_ke,linewidth = 5.0)
#        ax.yaxis.get_offset_text().set_fontsize(50)
#        plt.savefig(''.join(['plume_v3_stratified_keElement.png']))

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
        print("Particle Settling Velocity: {:2d}".format(particle_settlingVel))

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

#           if(particle_position_x > gridpoints_x):
#               particle_position_x = gridpoints_x
#               particle_x_velocity = 0
#           if(particle_position_x < 0):
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

        print(particle_position_vector)

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
        y = mp.mesh_generation(y_cluster,elements_y,y_start,y_end,order,1,'in')

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
            print("Iteration no.: %s" % k)

            # Reads data files.
            data,time,istep = rn.readnek(''.join([filename,'0.f',repr(k+1).zfill(5)]))

            # Compute timestep.
            dt = time - time_old
            time_old = time

            print("Time step:      %s" % dt)
            print("Actual time:    %s" % time)

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

        total_files = int(total_timesteps/jump)

        #file_counter = 1

        if (start_file == 1):
            range_vals = [x - (jump - 1) for x in np.array(range(1,total_files))*jump]
        else:
            range_vals = np.array(range(int(math.floor(start_file/jump)),total_files))*jump
            print("Make sure calculation of file numbers is correct")
            print(range_vals)

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
#               magVel = np.sqrt(np.square(verVel) + np.square(horVel))
                pressure = mesh[:,:,2]
            else:
#                temperature = mesh[:,:,5]
                horVel = mesh[:,:,2]
                verVel = mesh[:,:,3]
#               magVel = np.sqrt(np.square(verVel) + np.square(horVel))
                pressure = mesh[:,:,4]

            # Defines size of grid.
            x = mp.mesh_generation(x_cluster,elements_x,gridpoints_x,order,4,'out')
            y = mp.mesh_generation(y_cluster,elements_y,gridpoints_y,order,2,'out')

#           x_data = temperature[200,:]
#           name = '_tempLine_'
#           x1 = 0
#           x2 = 10
#           y1 = 0
#           y2 = 1
#           pt.myPlot(x,x_data,time,'Horizontal Position','Temperature',filename,name,k/jump,x1,x2,y1,y2)

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


def scatterPlume( filename, file_loc, order, dimension, start_file, jump, final_timestep, numPlots, elements_x, elements_y, elements_z ):

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

        order = order - 1

        final_file = int(final_timestep/jump)

        previous = 59
        ii = 0 + previous*2

        if (start_file == 1):
            range_vals = [x - (jump - 1) for x in np.array(range(1,final_file+1))*jump]

        else:
            range_vals = np.array(range(int(math.floor(float(start_file)/jump)),final_file+1))*jump
            print("Make sure calculation of file numbers is correct")
            print(range_vals)

        # Reading in mesh data.
        data,time,istep,header,elmap,u_i,v_i,w_i,t_i = rn.readnek(''.join([file_loc, \
                filename,'0.f00001']))

        [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

        # Defining x,y,z coordinates.
        x = mesh[:,0,0,0]
        y = mesh[0,:,0,1]
        z = mesh[0,0,:,2]

        x_start = x[0]
        x_end = x[-1]
        y_start = y[0]
        y_end = y[-1]
        z_start = z[0]
        z_end = z[-1]

        height = len(z)

        x_len = len(x)
        y_len = len(y)

        for t in range_vals:

            file_num = int((t-1)/jump + 1)

            # Outputs counter to terminal.
            if (start_file == 1):
                files_remaining = int(final_file - t/jump)
            else:
                files_remaining = int(final_file - (t-range_vals[0])/jump - start_file/jump)

            sys.stdout.write("\r")
            sys.stdout.write("Files remaining: {:2d}".format(files_remaining))
            sys.stdout.flush()

            # Reads data files.
            data,time,istep,header,elmap,u_i,v_i,w_i,t_i = \
                        rn.readnek(''.join([file_loc,filename,'0.f',repr(t).zfill(5)]))

            # Reshapes data onto uniform grid.
            if (dimension == 2):
                [ mesh ] = rn.reshapenek2D(data, elements_z, elements_x)

            elif (dimension == 3):

                [ mesh ] = rn.reshapenek3D(data, elements_x, elements_y, elements_z)

            temperature = mesh[:,:,:,t_i]

            x_loc3D = []
            y_loc3D = []
            z_loc3D = []

            # Running through heights
            for k in range(3,height):

                sum_total = 0
                x_tot = x_len - 1
                y_tot = y_len - 1

                tol = 0.1

#               x_pos = []
#               y_pos = []

                for i in range(0,x_tot):
                    for j in range(0,y_tot):
                        if temperature[i,j,k] > tol:

#                           x_pos.append(x[i])
#                           y_pos.append(y[j])

                            x_loc3D.append(x[i])
                            y_loc3D.append(y[j])
                            z_loc3D.append(z[k])

#               plt.figure(figsize=(20,20))
#               plt.scatter(x_pos,y_pos)#,outline)
#               plt.title(''.join(['Height = ',repr(z[k])]),fontsize=40)
#               plt.xlim(-3,3)
#               plt.ylim(-3,3)
#               plt.savefig(os.path.join(''.join(['height',repr(k).zfill(5),'.png']))\
#                                                                       ,bbox_inches='tight')
#               plt.close('all')


            pt.plot3DScatter(x_loc3D,y_loc3D,z_loc3D,x_start,x_end,y_start,y_end,\
                                                        z_start,z_end,order,ii,time,file_num,x,y,z)
            ii = ii + 2
            ii = ii % 360

        return
