#Link to path in which postProcess_lib is stored.
import sys
import os
sys.path.insert(1,'/home/ufaserv1_g/scdrw/Python')
import postProcess_lib as pp

# Chooses which postprocessing script to run.
# 0 - PseudoColourPlotting
# 1 - integrateDomain
# 2 - integratePlume

switch = 0

def pseudocolour():
        pp.PseudoColourPlotting( 'plume_v3_particles',
        7,      # Order 
        3,      # Dimension
        1,      # Start file
        10,     # Jump
        1954,   # Final timestep
        2,      # Number of plots
        80,     # Elements in x
        80,     # Elements in y
        10,     # Elements in z
        -5,     # x lower boundary
        5,      # x upper boundary
        0,      # y lower boundary
        10,     # y upper boundary
        0.25,   # Position of z slice
        0.98,   # Clustering in x
        0.98,   # Clustering in y
        0       # Particles (0 - no, 1 - yes)
        )
        return

def integrateDomain():
        pp.integrateDomain( 'plume_v2_3D',
        7,      # Order
        3,      # Dimension
        1,     # Jump
        1,      # Final timestep
        80,     # Elements in x
        40,     # Elements in y
        10,     # Elements in z
        -5,     # x lower boundary
        5,      # x upper boundary
        0,      # y lower boundary
        5,      # y upper boundary
        0,      # z lower boundary
        0.5,    # z upper boundary
        0.98,   # Clustering in x
        0.98,   # Clustering in y
        0       # Particles (0 - no, 1 - yes)
        )
        return

def integratePlume():
        pp.integratePlume( 'plume_v2_3D',
        7,      # Order
        3,      # Dimension
        1,     # Jump
        1,      # Final timestep
        80,     # Elements in x
        40,     # Elements in y
        10,     # Elements in z
        -5,     # x lower boundary
        5,      # x upper boundary
        0,      # y lower boundary
        5,      # y upper boundary
        0,      # z lower boundary
        0.5,    # z upper boundary
        0.98,   # Clustering in x
        0.98,   # Clustering in y
        0       # Particles (0 - no, 1 - yes)
        )
        return

def plotparticlesonly():
        pp.plotparticlesonly(
        7,      # Order 
        1,      # Start file
        50,     # Jump
        8140,   # Final particle file 
        80,     # Elements in x
        80,     # Elements in y
        10,     # Elements in z
        -5,     # x lower boundary
        5,      # x upper boundary
        0,      # y lower boundary
        10,     # y upper boundary
        0.25,   # Position of z slice
        0.98,   # Clustering in x
        0.98,   # Clustering in y
        )
        return

def choose_function(argument):
        switcher = {
                0: pseudocolour,
                1: integrateDomain,
                2: integratePlume,
                3: plotparticlesonly,
        }
        # Get the function from switcher dictionary
        func = switcher.get(argument)

        return func()

choose_function(switch)



#pp.meshPlot( 160, 160, 10, 10, 0.9, 0.92, 7 )

#pp.time_finder( 'plume_v9_meshInd', 1, 8000 )

# Combine images into video using the script make_video_from_image.sh
#os.system('./make_video_from_image.sh')

#pp.meshInd()

#pp.linePlot( 'plume_v2_fullDomain', 7, 1, 10, 1200, 80, 80, 10, 10, 0.9, 0.955555, 0 )
