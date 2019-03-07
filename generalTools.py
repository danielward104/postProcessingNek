import numpy as np

def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx


def mkdir_p(mypath):
        '''Creates a directory. equivalent to using mkdir -p on the command line'''

        from errno import EEXIST
        from os import makedirs,path
        
        try:
            makedirs(mypath)
        except OSError as exc: # Python >2.5
            if exc.errno == EEXIST and path.isdir(mypath):
                pass    
            else: raise
