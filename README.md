# postProcessingNek
Used to post-process data from a Nek5000 simulation.  The script postProcess_lib.py contains a number of functions used to convert the output of Nek5000 files to the correct format and then create line plots and pseudocolour plots.

The function readnek takes the Nek5000 output file (in binary) and reads it into a tensor.

The function reshapenek takes the output from readnek and uses the element map to map local elements back to their global positions on the grid.

The function geometricRatio is used to create a vector with all entries separated by geometric ratios.

The function PseudoColourPlotting plots data from a Nek5000 run.  Inputs are as follows:

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

The function integrateDomain plots line data from a Nek5000 run.  Inputs are as follows:

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

The functions myPcolour and myPlot are simply used to plot a pseudo-colour plot and a line plot respectively.
