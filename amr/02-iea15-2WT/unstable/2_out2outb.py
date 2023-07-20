import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.weio.fast_output_file import FASTOutputFile
from welib.tools.tictoc import Timer

filenames=[]
filenames+=['IEA-15-240-RWT-Monopile.T1.out']
filenames+=['IEA-15-240-RWT-Monopile.T2.T1.out']

for filename in filenames:
    print('Reading >>> ', filename)
    filenameOut = filename.replace('.out','.conv.outb')
    with Timer('read'):
        fo = FASTOutputFile(filename)
    print('Writting >>> ', filenameOut)
    with Timer('write'):
        fo.toOUTB(filename=filenameOut)

if __name__ == '__main__':
    pass
