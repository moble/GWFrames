#! /usr/bin/env python

"""
Convert old dat files to (NINJA-style) h5 format.

"""

if __name__ == "__main__" :
    import os
    import os.path
    import sys
    import re
    import glob
    import argparse
    import numpy
    import h5py

    # Set up and run the parser
    parser = argparse.ArgumentParser(description = __doc__)
    # parser.add_argument('--', action='store_true',
    #                     help='')
    parser.add_argument('--input_directory', default='./',
                        help='Directory to search for dat files.')
    parser.add_argument('--output_directory', default='./',
                        help='Directory in which to place output files.')
    parser.add_argument('--ADMMass',
                        help='ADM Mass to place in the output files.')
    args = vars(parser.parse_args(sys.argv[1:]))

    LapseFiles = glob.glob(args['input_directory']+'/LapseSurfaceIntegral_R*.dat')
    hFiles = glob.glob(args['input_directory']+'/rh_R*m.dat')
    Psi4Files = glob.glob(args['input_directory']+'/rPsi4_R*m_U8+.dat')
    AreaFiles = glob.glob(args['input_directory']+'/SurfaceArea_R*.dat')

    h = h5py.File(args['output_directory']+'/rh_FiniteRadii_CodeUnits.h5', 'w')
    Psi4 = h5py.File(args['output_directory']+'/rPsi4_FiniteRadii_CodeUnits.h5', 'w')

    for LapseFile,hFile,Psi4File,AreaFile in zip(LapseFiles,hFiles,Psi4Files,AreaFiles) :
        LapseData = numpy.loadtxt(LapseFile)
        hData = numpy.loadtxt(hFile)
        Psi4Data = numpy.loadtxt(Psi4File)
        AreaData = numpy.loadtxt(AreaFile)
        R = re.search('LapseSurfaceIntegral_(.*)m.dat', LapseFile).group(1)
        print(R)
        LMData = []
        with open(hFile, 'r') as f :
            line = f.readline()
            while(line.startswith('#')) :
                match = re.search('.*\((.*?),(.*?)\)\(R=.*', line)
                if(match) :
                    LMData.append([int(ellm) for ellm in match.groups()])
                line = f.readline()
        RData_h = h.create_group(R+".dir")
        RData_Psi4 = Psi4.create_group(R+".dir")
        ArealRadius = numpy.sqrt(AreaData[:,1]/(4*numpy.pi))
        AverageLapse = LapseData[:,1]/AreaData[:,1]
        RData_h.create_dataset('ArealRadius.dat', data=numpy.array([AreaData[:,0], ArealRadius]).transpose(),
                               compression="gzip", shuffle=True)
        RData_h.create_dataset('AverageLapse.dat', data=numpy.array([LapseData[:,0], AverageLapse]).transpose(),
                               compression="gzip", shuffle=True)
        RData_h.create_dataset('CoordRadius.dat', data=[[0., float(R[1:])]],
                               compression="gzip", shuffle=True)
        RData_h.create_dataset('InitialAdmEnergy.dat', data=[[0., float(args['ADMMass'])]],
                               compression="gzip", shuffle=True)
        RData_Psi4.create_dataset('ArealRadius.dat', data=numpy.array([AreaData[:,0], ArealRadius]).transpose(),
                                  compression="gzip", shuffle=True)
        RData_Psi4.create_dataset('AverageLapse.dat', data=numpy.array([LapseData[:,0], AverageLapse]).transpose(),
                                  compression="gzip", shuffle=True)
        RData_Psi4.create_dataset('CoordRadius.dat', data=[[0., float(R[1:])]],
                                  compression="gzip", shuffle=True)
        RData_Psi4.create_dataset('InitialAdmEnergy.dat', data=[[0., float(args['ADMMass'])]],
                                  compression="gzip", shuffle=True)
        for i,(ell,m) in enumerate(LMData) :
            if(i%2) :
                Ylm = RData_h.create_dataset("Y_l{0}_m{1}.dat".format(ell, m), data=numpy.array([hData[:,0], hData[:,i], hData[:,i+1]]).transpose(),
                                             compression="gzip", shuffle=True)
                Ylm.attrs['ell'] = ell
                Ylm.attrs['m'] = m
                Ylm = RData_Psi4.create_dataset("Y_l{0}_m{1}.dat".format(ell, m), data=numpy.array([Psi4Data[:,0], Psi4Data[:,i], Psi4Data[:,i+1]]).transpose(),
                                                compression="gzip", shuffle=True)
                Ylm.attrs['ell'] = ell
                Ylm.attrs['m'] = m

    h.close()
    Psi4.close()
