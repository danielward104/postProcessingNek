import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np

# Plots pseudocolour.
def myPcolour(x,y,data,time,x_label,y_label,filename,name,file_counter,**kwargs):

        domain_x = x[0,-1] - x[0,0]
        domain_y = y[-1,0] - y[0,0]

        if (domain_y - domain_x > 0):
            ratio = domain_x/domain_y
            domain_y = 25
            domain_x = ratio*25
        else:
            ratio = domain_y/domain_x
            domain_x = 25
            domain_y = ratio*25

        plt.figure(figsize=(domain_x, domain_y)) # Increases resolution.
	plt.title(''.join([name,', time = %2d'%round(time,2)]),fontsize=40)
#        plt.title('time = %2d'%(time),fontsize=40)
        plt.xlabel(x_label,fontsize=40)
        plt.ylabel(y_label,fontsize=40)
        plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)
	plt.pcolormesh(x,y,data,**kwargs)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize = 30)  # vertically oriented colorbar
        plt.savefig(''.join([filename,'_',name,'_',repr(file_counter).zfill(5),'.png']),bbox_inches='tight')
#       plt.savefig('temp.png')

        plt.close('all')

        return

# Plots pseudocolour, including particle behaviour.
def particlePcolour(x,y,data,time,x_label,y_label,filename,name,file_counter,x_ppos,y_ppos,**kwargs):

        particlesOnly = 0

        if(particlesOnly < 1):
            plt.figure(figsize=(50, 25)) # Increases resolution.
            plt.title(''.join([name,', time = %2d'%round(time,2)]),fontsize=40)
            plt.xlabel(x_label,fontsize=40)
            plt.ylabel(y_label,fontsize=40)
            plt.xticks(fontsize = 30)
            plt.yticks(fontsize = 30)
            plt.pcolormesh(x,y,data,**kwargs)
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize = 30)  # vertically oriented colorbar

        #plt.scatter(x_ppos,y_ppos,marker='.',color='black',s=50)

        #if(particlesOnly > 0):
        #plt.axis([30,50,0,40])

        plt.savefig(''.join([filename,'_',name,'_',repr(file_counter).zfill(5),'_particle.png']))

        plt.close('all')

        return

# Plots pseudocolour, with velocity quiver plot.
def myPcolourQuiver(x,y,data,quiver_x,quiver_y,time,x_label,y_label,filename,name,file_counter,**kwargs):

	domain_x = x[0,-1] - x[0,0]
	domain_y = y[-1,0] - y[0,0]

	if (domain_y - domain_x > 0):
	    ratio = domain_x/domain_y
	    domain_y = 25
	    domain_x = ratio*25
	else:
	    ratio = domain_y/domain_x
	    domain_x = 25
	    domain_y = ratio*25

        plt.figure(figsize=(int(domain_x) + 10, int(domain_y))) # Increases resolution.
        plt.title(''.join([name,', time = %2d'%round(time,2)]),fontsize=40)
#        plt.title('time = %2d'%(time),fontsize=40)
        plt.xlabel(x_label,fontsize=40)
        plt.ylabel(y_label,fontsize=40)
        plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)
        plt.pcolormesh(x,y,data,**kwargs)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize = 30)  # vertically oriented colorbar

	x_length = len(np.transpose(x))
	y_length = len(y)

	x_quiver = np.zeros((y_length,x_length))

	for i in range(0,y_length):
	    x_quiver[i,:] = x

        y_quiver = np.zeros((y_length,x_length))
        for i in range(0,x_length):
            y_quiver[:,i] = np.transpose(y)

	scale = 7

	magVel = np.sqrt(np.square(quiver_x[::scale,::scale]) + np.square(quiver_y[::scale,::scale]))
	plot_u = quiver_x[::scale,::scale]/magVel
	plot_v = quiver_y[::scale,::scale]/magVel
	where_are_NaNs = np.isnan(plot_u)
	plot_u[where_are_NaNs] = 0
	where_are_NaNs = np.isnan(plot_v)
	plot_v[where_are_NaNs] = 0

        plt.quiver(x[::scale,::scale],y[::scale,::scale],plot_u,plot_v, magVel,cmap='RdBu_r',scale=50,width=0.001)

        plt.savefig(''.join([filename,'_',name,'_',repr(file_counter).zfill(5),'.png']),bbox_inches='tight')
#       plt.savefig('temp.png')

        plt.close('all')

        return

# Plots line plots.
def myPlot(x,y,time,x_label,y_label,filename,name,file_counter,x1,x2,y1,y2,orientation):

	if(orientation == 'long'):
	    plt.figure(figsize=(25, 15)) # Increases resolution.
	    yplot = np.zeros(len(x))
            plt.plot(x,yplot,color='black',linewidth=0.5)
	elif(orientation == 'thin'):
	    plt.figure(figsize=(15, 25)) # Increases resolution.
	    xplot = np.zeros(len(y))
            plt.plot(xplot,y,color='black',linewidth=0.5)

	xplot = np.zeros(len(y))
        plt.plot(xplot,y,color='black',linewidth=0.5)

	plt.title('time = %d'%round(time,3),fontsize=40)
        plt.xlabel(x_label,fontsize=40)
        plt.ylabel(y_label,fontsize=40)
	plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)
	plt.xlim(x1,x2)
	plt.ylim(y1,y2)
        plt.plot(x,y)

        plt.savefig(''.join([filename,name,repr(file_counter).zfill(5),'.png']))

        plt.close('all')

        return

