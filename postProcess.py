#Link to path in which postProcess_lib is stored.
import sys
import os
sys.path.insert(1,'/home/cserv1_a/soc_pg/scdrw/Documents/nbudocuments/PhD/SimNumerics/Python/postProcessingLib')
import postProcess_lib as pp

pp.PseudoColourPlotting( 'plume_v2_3D', 
7,	# Order 
3, 	# Dimension
2, 	# Start file
1, 	# Jump
2, 	# Final timestep
1, 	# Number of plots
80, 	# Elements in x
40, 	# Elements in y
10, 	# Elements in z
-5, 	# x lower boundary
5, 	# x upper boundary
0, 	# y lower boundary
5,	# y upper boundary
0.25, 	# Position of z slice
0.98, 	# Clustering in x
0.98, 	# Clustering in y
0 	# Particles (0 - no, 1 - yes)
)

#def PseudoColourPlotting( filename, order, dimension, start_file, jump, total_timesteps, numPlots, elements_x, elements_y, elements_z, x_start, x_end, y_start, y_end, z_slice, x_cluster, y_cluster, gridType, particles ):

#pp.meshPlot( 160, 160, 10, 10, 0.9, 0.92, 7 )

#pp.time_finder( 'plume_v9_meshInd', 1, 8000 )

#pp.integrateDomain( 'plume_v3_stratified', 200, 15000, 40, 80, 50, 100, 0.8, 1.1, 0 )

# Combine images into video using the script make_video_from_image.sh
#os.system('./make_video_from_image.sh')

#pp.meshInd()

#pp.linePlot( 'plume_v2_fullDomain', 7, 1, 10, 1200, 80, 80, 10, 10, 0.9, 0.955555, 0 )
