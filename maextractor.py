import sys
import os
import traceback
import csv
import pandas
import numpy
import matplotlib.pyplot as plt
from fisher_py import RawFile
from fisher_py.data.business import TraceType
from scipy.integrate import trapezoid
from pathlib import Path

#List within which the MAs are held
MAlist=[]
#inputread = input("What file do you want to read in? ")

# Reading in and validating the requested file
def extractor(readIn):
    raw_file = RawFile(readIn)

    # Setting the parameters.
    target_mass = 137.0464
    mass_tolerance_ppm = 10
    rt, i = raw_file.get_chromatogram(target_mass, mass_tolerance_ppm, TraceType.MassRange)
    #mz, i2, charges, real_rt = raw_file.get_scan_ms1(1)

    # plt.figure()
    # graph = plt.plot(rt, i)
    # plt.savefig('myplot1.png')

    # Gets the peak so that we know which one to integrate
    peak_index = numpy.argmax(i)
    peak_rt= rt[peak_index]
    baseline = 0.0

    # Calculates the left boundary by finding where you get back to baseline before the peak
    ####### This needs to be refined using the savgol filter first, because that's what xcalibur does to get rid of noise
    intensity = numpy.argmax(i)
    leftbound = peak_index
    while (intensity > baseline and intensity-i[leftbound+1] < 50000):
        leftbound -= 1
        intensity = i[leftbound]

    # Calculates the right boundary by finding where you get back to baseline before the peak
    intensity = numpy.argmax(i)
    rightbound = peak_index
    while (intensity > baseline and intensity - i[rightbound-1] < 50000): # Ask what sensitivity they want for this
        rightbound += 1
        intensity = i[rightbound]


    # Makes a slice of the retention time and intensity lists so that it only involves the peak's curve
    peakslicex = rt[leftbound:rightbound]
    peakslicey =i[leftbound:rightbound]

    # Calculates the MA value (area under the peak's curve)
    MAval = trapezoid(peakslicey, peakslicex)
    MAlist.append(MAval)

def filepather():
    # dirinput = input("What's the name of the target directory? ")
    outputwrite = input("What file do you want to output? ")

    path = Path.cwd()
    #q = path / dirinput

    pathlist = []
    samplenames = []

    for child in path.iterdir():
        if child.name.endswith(".raw"):
            pathlist.append(child)

    for p in pathlist:
        rightpath = Path(p).name
        samplenames.append(rightpath.removesuffix('.raw'))
        extractor(rightpath)

        # Saves the MAlist to a csv file that you specify
    with open(outputwrite, 'a', newline='') as csvfile:
        dawriter = csv.writer(csvfile)
        dawriter.writerow(samplenames)
        dawriter.writerow(MAlist)

filepather()

