# Link to path in which postProcess_lib is stored.
import sys
import os
sys.path.insert(1,'/home/cserv1_a/soc_pg/scdrw/Documents/nbudocuments/PhD/SimNumerics/Python/postProcessingLib')
import postProcess_lib as pp

pp.PseudoColourPlotting( 'plume_v2_fullDomain', 7, 1, 20, 1156, 4, 60, 60, 10, 10, 0.9, 0.9, 1, 0 )

#pp.meshPlot( 60, 60, 10, 10, 0.9, 0.92, 7 )

#pp.time_finder( 'plume_v9_meshInd', 1, 8000 )

#pp.integrateDomain( 'plume_v9_meshInd', 200, 15000, 40, 80, 50, 100, 0.8, 1.1, 0 )

# Combine images into video using the script make_video_from_image.sh
#os.system('./make_video_from_image.sh')

#pp.meshInd()
