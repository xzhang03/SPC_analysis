from sdtfile import SdtFile
import numpy as np

def readfile(path):
    sdt = SdtFile(path)
    reshaped=np.reshape(sdt.data[0],(512,-1,256))[:,:1250,:]
    return reshaped