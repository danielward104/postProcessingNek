import matplotlib.pyplot as plt

# Plots pseudocolour.
def myPcolour(x,y,data,time,x_label,y_label,x_range,y_range,filename,name,file_counter,**kwargs):

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
#       plt.savefig('temp.png')

        plt.close('all')

        return

# Plots pseudocolour, including particle behaviour.
def particlePcolour(x,y,data,time,x_label,y_label,x_range,y_range,filename,name,file_counter,x_ppos,y_ppos,**kwargs):

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

        #plt.scatter(x_ppos,y_ppos,marker='.',color='black',s=50)

        #if(particlesOnly > 0):
        #plt.axis([30,50,0,40])

        plt.savefig(''.join([filename,name,repr(file_counter).zfill(5),'_particle.png']))

        plt.close('all')

        return

# Redundant?
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

