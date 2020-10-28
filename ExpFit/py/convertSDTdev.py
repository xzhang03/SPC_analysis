import sdtfile
from scipy import io
import numpy as np
import os
import sys
import glob

rootdir = "H:/2p/stephen"
sdtFiles = glob.glob(rootdir+'/*.sdt',recursive=True)
print(os.path.isdir('H:/2p'))