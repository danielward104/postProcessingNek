import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np

# Plots pseudocolour.
def myPcolour(x,y,data,time,xmin,xmax,ymin,ymax,x_label,y_label,filename,name,file_counter,**kwargs):

        domain_x = xmax - xmin
        domain_y = ymax - ymin

        if (domain_y - domain_x > 0):
            ratio = domain_x/domain_y
            domain_y = 25
            domain_x = ratio*25
        else:
            ratio = domain_y/domain_x
            domain_x = 25
            domain_y = ratio*25

        plt.figure(figsize=(domain_x, domain_y)) # Increases resolution.
        plt.title(''.join([name,', time = %5.3f'%(time)]),fontsize=40)
        axes = plt.gca()
        axes.set_xlim([xmin,xmax])
        axes.set_ylim([ymin,ymax])
        plt.xlabel(x_label,fontsize=40)
        plt.ylabel(y_label,fontsize=40)
        plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)

        # Pseudocolour
        #plt.pcolormesh(x,y,data,**kwargs)

        # Contour plot
        #plt.contour(x,y,data,100,**kwargs)

        # Filled contour plot
        plt.contourf(x,y,data,100,**kwargs)
        cbar = plt.colorbar()
        #plt.contour(x,y,data,100,colors='k')

        cbar.ax.tick_params(labelsize = 30)  # vertically oriented colorbar
        plt.savefig(''.join([filename,'_',name,'_',repr(file_counter).zfill(5),'.png']), \
            bbox_inches='tight')
#       plt.savefig('temp.png')

        plt.close('all')

        return

# Plots pseudocolour, including particle behaviour.
def particlePcolour(x,y,data,time,x_label,y_label,filename,name,file_counter,x_ppos,y_ppos,**kwargs):

#       domain_x = x[0,-1] - x[0,0]
#        domain_y = y[-1,0] - y[0,0]
#
#        if (domain_y - domain_x > 0):
#            ratio = domain_x/domain_y
#            domain_y = 25
#            domain_x = ratio*25
#        else:
#            ratio = domain_y/domain_x
#            domain_x = 25
#            domain_y = ratio*25

        domain_x = 20
        domain_y = 10

        plt.figure(figsize=(domain_x, domain_y)) # Increases resolution.
        plt.title(''.join([name,', time = %2d'%round(time,2)]),fontsize=40)
        plt.xlabel(x_label,fontsize=40)
        plt.ylabel(y_label,fontsize=40)
        plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)
        plt.pcolormesh(x,y,data,**kwargs)

        axes = plt.gca()
        axes.set_xlim([min(x),max(x)])
        axes.set_ylim([min(y),max(y)])

        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize = 30)  # vertically oriented colorbar

        plt.scatter(x_ppos,y_ppos,marker='.',color='black',s=5)

        plt.savefig(''.join([filename,'_',name,'_',repr(file_counter).zfill(5),'_particle.png']), \
            bbox_inches='tight')

        plt.close('all')

        return


# Plots particles only.
def particleOnlyPlot(x,y,x_label,y_label,file_counter,x_ppos,y_ppos,):

        domain_x = 20
        domain_y = 10

        plt.figure(figsize=(domain_x, domain_y)) # Increases resolution, decreases plotting speed.
        plt.title(''.join(['File number = %2d'%file_counter]),fontsize=40)
        plt.xlabel(x_label,fontsize=40)
        plt.ylabel(y_label,fontsize=40)
        plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)

        axes = plt.gca()
        axes.set_xlim([min(x),max(x)])
        axes.set_ylim([min(y),max(y)])

        plt.scatter(x_ppos,y_ppos,marker='.',color='black',s=1)

        plt.savefig(''.join(['particle_',repr(file_counter).zfill(5),'.png']), \
            bbox_inches='tight')

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

def meshPlot(x,y):

        domain_x = x[-1] - x[0]
        domain_y = y[-1] - y[0]

        if (domain_y - domain_x > 0):
            ratio = domain_x/domain_y
            domain_y = 40
            domain_x = ratio*40
            
        else:
            ratio = domain_y/domain_x
            domain_x = 40
            domain_y = ratio*40

        plt.figure(figsize=(domain_x, domain_y)) # Increases resolution.
        plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)

        for i in range(0,len(x)):
            xplot = np.zeros(len(y))
            xplot = [q + x[i] for q in xplot]
            plt.plot(xplot,y,color='black',linewidth=1)
        
        for j in range(0,len(y)):
            yplot = np.zeros(len(x))
            yplot = [p + y[j] for p in yplot]
            plt.plot(x,yplot,color='black',linewidth=1)

        plt.xlim(x[0],x[-1])
        plt.ylim(y[0],y[-1])

        plt.savefig('mesh.jpg',bbox_inches='tight')

        return





