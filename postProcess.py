#Link to path in which postProcess_lib is stored.
import sys
import os
sys.path.insert(1,'/home/cserv1_a/soc_pg/scdrw/Documents/nbudocuments/PhD/SimNumerics/Python/postProcessingLib')
import postProcess_lib as pp

# Chooses which postprocessing script to run.
# 0 - PseudoColourPlotting
# 1 - integrateDomain
# 2 - integratePlume

switch = 2 

def pseudocolour():
	pp.PseudoColourPlotting( 'plume_v9_axi', 
	7,	# Order 
	3, 	# Dimension
	1, 	# Start file
	1, 	# Jump
	1, 	# Final timestep
	2, 	# Number of plots
	8, 	# Elements in x
	40, 	# Elements in y
	8, 	# Elements in z
	0.0, 	# Position of z slice
	0 	# Particles (0 - no, 1 - yes)
	)
	return	

def integrateDomain():
	pp.integrateDomain( 'plume_v4_LES',         
	7,     	# Order
	3,      # Dimension
	10,	# Jump
	100,	# Final timestep
	12,     # Elements in x
	6,     # Elements in y
	3,     # Elements in z
	-2,     # x lower boundary
	2,      # x upper boundary
	0,      # y lower boundary
	2,      # y upper boundary
	0,      # z lower boundary
	1.0,    # z upper boundary
	0.9998,   # Clustering in x
	0.9998,   # Clustering in y
	0       # Particles (0 - no, 1 - yes)
	)
	return

def integratePlume():
        pp.integratePlume( 'plume_v9_axi',
        7,      # Order
        0.003,  # Tolerance
        3,      # Dimension
        50,     # Jump
        2466,	# Final timestep
        8 ,     # Elements in x
        20,     # Elements in y
        8,	# Elements in z
        0.0,    # Position of z slice
        0.45,    # Umbrella cutoff
        )
        return

def average_field():
        pp.average_field('plume_v9_axi',
        7,      # Order 
        3,      # Dimension
        3000,      # Start file
        5,      # Jump
        7000,      # Final timestep
        8,      # Elements in x
        30,     # Elements in y
        8      # Elements in z
        )
        return

def TKE():
        pp.TKE('plume_v9_axi',
        7,      # Order 
        3,      # Dimension
        3000,      # Start file
        5,      # Jump
        6113,      # Final timestep
        8,      # Elements in x
        30,     # Elements in y
        8      # Elements in z
        )
        return
	
def choose_function(argument):
	switcher = {
		0: pseudocolour,
		1: integrateDomain,
	        2: integratePlume,
                3: average_field,
                4: TKE,
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
