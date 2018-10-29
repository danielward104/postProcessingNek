import numpy as np
import math
from array import array
import sys

# Module to read the Nek5000 data files of the form 0.f##### and convert them to a spectral grid.

def readnek( fname ):

        #--------------------------------------------
        # Open the file
        #--------------------------------------------
    arr = array('L')
    with open(fname,'rb') as f:

        #--------------------------------------------
        # Read header of tile
        #--------------------------------------------
        header = f.readline(132)
        etag = f.read(4)
        # Precision of simulation. 
        wdsz = float(header[5])

        # Element sizes
        lr1 = [float(header[7:9].decode("utf-8").strip()), \
            float(header[10:13].decode("utf-8").strip()),float(header[13:16].decode("utf-8").strip())]

        # Compute total number of points per element
        npel = 1
        for i in range(0,len(lr1)):
            npel = npel*lr1[i]
        npel = int(npel)

        # Compute number of active dimensions
        if (lr1[2] > 1):
            add_dims = 1
        else:
            add_dims = 0
        ndim = 2 + add_dims

        # Number of elements
        nel = float(header[16:26].decode("utf-8").strip())

        # Number of elements in the file
        nelf = int(header[27:37].decode("utf-8").strip())

        # Time
        time = float(header[38:58].decode("utf-8").strip())

        # Iteration number
        istep = float(header[59:68].decode("utf-8").strip())

        # Get fields [XUPT]
        fields = header[83:].decode("utf-8")
        var = np.zeros(5)
        if (int('X' in fields) == 1):
            var[0] = ndim
        if (int('U' in fields) == 1):
            var[1] = ndim
        if (int('P' in fields) == 1):
            var[2] = 1
        if (int('T' in fields) == 1):
            var[3] = 1

        # Saves indices of velocity and temperature.
        u_i = (ndim + 1) - 1
        v_i = (ndim + 2) - 1
        w_i = (ndim + 3) - 1
        t_i = (2*ndim + 2) - 1

        # This will probably need resorting for 2D postprocessing, re. w_i and ndim = 2.

        # Gets number of fields
        nfields = 0
        for i in range(0,len(var)):
            nfields = nfields + var[i]
        nfields = int(nfields)

        # Read element map
        #map_proxy = f.readlines()[1:]
        elmap = np.fromfile(f,dtype='int32',count=nelf)

        #--------------------------------------------
        # Read data
        #--------------------------------------------
        data = np.zeros((nelf,npel,nfields))
        for ivar in range(1,6):
            if (ivar == 1):
                idim0 = 0
            else:
                idim0 = 0
                for i in range(0,ivar-1):
                    idim0 = idim0 + var[i]
            for iel in elmap:
                iter_range = [x+idim0 for x in range(1,int(var[ivar-1])+1)]
                iter_range = [int(x) for x in iter_range]
                for idim in iter_range:
                    data[iel-1,:,idim-1] = np.fromfile(f,dtype='float32',count=npel)

    return [data,time,istep,header,elmap,u_i,v_i,w_i,t_i]


def reshapenek3D( data, nelx, nely ,nelz ):
        nel,N3,nfld = data.shape
        N = round(math.pow(float(N3),1.0/3.0),0)

        if (nel != nelx*nely*nelz):
            print('Error in reshapenek: nel != nelx*nely*nelz.')
            sys.exit()

        #--------------------------------------------
        # Reshape data
        #--------------------------------------------

        mesh = np.zeros((int((N-1)*nelx+1),int((N-1)*nely+1),int((N-1)*nelz+1),nfld))

        for ifld in range(0,nfld):

            for iel in range(0,nel):

                ielz = math.floor(iel/(nelx*nely)) + 1
                iely = (math.floor(iel/nelx) % nely) + 1
                ielx = (iel % nelx) + 1

                ii = [x+(N-1)*(ielx-1) for x in range(0,int(N))]
                ii = [int(x) for x in ii]
                jj = [y+(N-1)*(iely-1) for y in range(0,int(N))]
                jj = [int(y) for y in jj]
                kk = [z+(N-1)*(ielz-1) for z in range(0,int(N))]
                kk = [int(z) for z in kk]

                mesh[ii[0]:(ii[7]+1),jj[0]:(jj[7]+1),kk[0]:(kk[7]+1),ifld] = \
                        np.reshape(data[iel,:,ifld], (8,8,8)).transpose()

        return [ mesh ]


def reshapenek3D_single( data, nelx, nely ,nelz ):
        nel,N3,nfld = data.shape
        N = round(math.pow(float(N3),1.0/3.0),0)

        if (nel != nelx*nely*nelz):
            print('Error in reshapenek: nel != nelx*nely*nelz.')
            sys.exit()

        #--------------------------------------------
        # Reshape data
        #--------------------------------------------

        ifld = 4
        z_slice = 4

        mesh = np.zeros((int((N-1)*nelx+1),int((N-1)*nely+1),int((N-1)*nelz+1)))

        for iel in range(0,nel):

            ielz = math.floor(iel/(nelx*nely)) + 1
            iely = (math.floor(iel/nelx) % nely) + 1
            ielx = (iel % nelx) + 1

            ii = [x+(N-1)*(ielx-1) for x in range(0,int(N))]
            ii = [int(x) for x in ii]
            jj = [y+(N-1)*(iely-1) for y in range(0,int(N))]
            jj = [int(y) for y in jj]
            kk = [z+(N-1)*(ielz-1) for z in range(0,int(N))]
            kk = [int(z) for z in kk]

            mesh[ii[0]:(ii[7]+1),jj[0]:(jj[7]+1),kk[0]:(kk[7]+1)] = \
                np.reshape(data[iel,:,ifld], (8,8,8)).transpose()

            mesh = mesh[:,:,z_slice]

        return [ mesh ]



def reshapenek2D( data, nelx, nely ):
        nel,N2,nfld = data.shape
        N = math.sqrt(N2)
        if (nel != nelx*nely):
            print('Error in reshapenek: nel != nelx*nely.')
            sys.exit()

        #--------------------------------------------
        # Reshape data
        #--------------------------------------------

        mesh = np.zeros((int((N-1)*nelx+1),int((N-1)*nely+1),nfld))

        for ifld in range(0,nfld):
            # Check this mesh isn't transposed!!!
#           mesh = np.zeros((int((N-1)*nelx+1),int((N-1)*nely+1)))

            for iel in range(0,nel):
                ielx = math.floor(iel/nely) + 1
                iely = (iel % nely) + 1

                ii = [x+(N-1)*(ielx-1) for x in range(0,int(N))]
                ii = [int(x) for x in ii]
                jj = [x+(N-1)*(iely-1) for x in range(0,int(N))]
                jj = [int(x) for x in jj]

                mesh[ii[0]:(ii[7]+1),jj[0]:(jj[7]+1),ifld] = np.reshape(data[iel,:,ifld], (8,8))

        return [ mesh ]

