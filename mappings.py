import numpy as np
import math

#def geometric_ratio( r, n, Sn ):
        # Calculating the axis using a geometric ratio.  The variable r is the geometric ratio to be used, n is the number of elements and Sn is the length of the axis, clusterEdge chooses which side to cluster the gridpoints on.

        # Compute first step size.
#        a = (r - 1)*Sn/(r**n - 1)

#        x = None
#        geosum = 0
#        geosum_last = 0

#        for i in range(0,n):

#            geosum = geosum + a*r**(i)

            # Computing spacing for each element
#            small_vector = np.linspace(geosum_last,geosum,8)

            # Concatenating vectors.
#            if (i < n-1):
#                small_vector = small_vector[0:7]
#            if (i == 0):
#                x = small_vector
#            else:
#                x = np.concatenate([x, small_vector])
#            geosum_last = geosum

#        return [ x ]

# Computes positions of nodes for SEM polynomials.
def lglnodes( N ):

        N1 = N + 1

        x = np.zeros((N+1,1))

        for i in range(0,N+1):
            x[i] = np.cos(math.pi*i/N)

        x = np.transpose(x)

        P = np.zeros((N1,N1))

        eps = 2.2204e-16
        check = 2

        while (np.amax(check) > eps):

            x_old = x

            P[:,0] = 1
            P[:,1] = x

            for k in range(2,N+1):
                 P[:,k] = ( (2*k - 1)*np.multiply(x,P[:,k-1]) - (k-1)*P[:,k-2] )/(k)

            x = x_old - ( np.multiply(x,P[:,N1-1]) - P[:,N-1] )/( N1*P[:,N1-1] )

            check = [abs(i) for i in x - x_old]

        return x

# Maps nodes calculated on a local element (in lglnodes) to a global element.
def mapping( y2, y1, order ):

        roots = lglnodes(order)

        roots_new = np.ones(order-1)*(y1 + y2)/2 + roots[0,1:order]*(y2 - y1)/2

        return roots_new

# Calculating the axis using a geometric ratio and Legendre polynomials.  The variable r is the geometric ratio to be used, n is the number of elements and Sn is the length of the axis, clusterEdge chooses which side to cluster the gridpoints on.
#def mesh_generation( r, Nel, Sn, order ):

#        mesh = np.zeros((order*Nel+1,1))
#        mesh[order] = (r - 1)*Sn/(r**Nel - 1)

#        x = mapping(mesh[0],mesh[order],order)

#        for i in range(1,order):
#            mesh[i] = x[i-1]

#        for j in range(2,Nel+1):
#            i = order*j
#            mesh[i] = mesh[i-order] + mesh[order]*r**(j-1)

#            x = mapping(mesh[i-order],mesh[i],order)
#            for k in range(1,order):
#                mesh[i-order+k] = x[k-1]

#        return mesh


# Generates nodes for the SEM.
def SEM_nodes( domain, Nel, N, cluster_dir, cluster_val ):

	# Input parameters:
	#domain = 10             # Size of the domain.
	#Nel = 52                # Number of elements.
	#N = 2                   # Number of clusters.
	#cluster_dir = 'out'     # Direction in which to cluster (out: clusters near edges).
	#cluster_val = 0.9       # Amount to cluster elements (geometric ratio).

	domain = float(domain)
	
	if (cluster_dir == 'out'):
	    cluster = 1/cluster_val
	else:
	    cluster = cluster_val

	for ix in range(0,N):

	    x = [0]*(Nel/N + 1)

	    if (ix == 0):
	        ratio = cluster
	    elif (ix == 1):
	        ratio = 1/cluster
	    elif (ix == 2):
	        ratio = cluster
	    else:
	        ratio = 1/cluster

	    nelx = Nel/N
	    x0 = ix*domain/N
	    x1 = (ix+1)*domain/N
	
	    dx = 1.0
	    x[0] = x0
	    for e in range(1,nelx+1):
	        x[e] = x[e-1] + dx
	        dx = ratio*dx

	    xlength = x[nelx] - x[0]
	    scale = (x1-x0)/xlength

	    for e in range(1,nelx+1):
	        x[e] = x0 + (x[e]-x0)*scale

	    if (ix > 0):
		x = x[1:]
	        x_all = x_all + x
	    else:
		x_all = x
		
	    for i in range(0,int(nelx/6)+1):
	        fir = 6*i
	        sec = 6*(i + 1)
	        if  (ix > 0):
	            fir = fir + 1
	            sec = sec + 1
	        if (sec > nelx):
	            sec = nelx + 1
	
	return x_all	

# Calculating the axis using a geometric ratio and Legendre polynomials.  The variable r is the geometric ratio to be used, n is the number of elements and Sn is the length of the axis, clusterEdge chooses which side to cluster the gridpoints on.
def mesh_generation( r, Nel, Sn, order, N, cluster_dir ):

	elements = SEM_nodes( Sn, Nel, N, cluster_dir, r )

        mesh = np.zeros((order*Nel+1,1))
	for i in range(0,Nel+1):
	    j = i*order
            mesh[j] = elements[i]
	
        x = mapping(mesh[0],mesh[order],order)

        for i in range(1,order):
            mesh[i] = x[i-1]

        for j in range(2,Nel+1):
            i = order*j
            x = mapping(mesh[i-order],mesh[i],order)
            for k in range(1,order):
                mesh[i-order+k] = x[k-1]

        return mesh
	
