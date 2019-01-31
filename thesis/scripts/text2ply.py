#!/usr/bin/env python3
from sys import argv

if (__name__=="__main__"):
    if (len(argv)<3):
        print("Usage: text2ply.py 'file with ascii tabulardata' outfile.ply")
        print("This program converts ascii data to a .ply file.")
        exit()
    
    import numpy as np
    import glob, os

    from third_party.python_plyfile.plyfile import PlyElement, PlyData

    ifile = argv[1];
    ofile = argv[2];

    if (not os.path.isfile(ifile)):
        print("Argument is not a file nor directory")
        exit()

    if (os.path.isdir(ifile)):
        ipath = ifile
        # load all data.txt files and write them to one output file
        #os.chdir(ipath)

        print("Not implemented yet");
        exit()

        i = 1;
        for ifile in glob.glob(ipath + "/*.data.txt"):
            idata = np.loadtxt(ifile)
            idata_d = [idata[:,:2]]
            if (i == 1):
                data = idata_d
                labels = i
            else:
                data = np.concatenate((data, idata_d), axis=0)
                labels = np.append(labels, i)
            i += 1
       
        print("Print to", argv[2])
        print("Total number of data samples:", np.size(labels))
        print("Total number of data points:", np.size(data))
    else:
        print("Load", ifile);
        data = np.loadtxt(ifile, usecols=(0,1,2))
        x = data[:,0];
        y = data[:,1];
        z = data[:,2];
        dp = np.array([x,y,z]);
        l0 = (x[0],y[0],z[0]);
        t = [tuple(row) for row in data]

        vertex = np.array(t, dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4')])

        el = PlyElement.describe(vertex, 'vertex')
        print("Write", ofile);

        with open(ofile, mode='wb') as f:
            PlyData([el]).write(f)

