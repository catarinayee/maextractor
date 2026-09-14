import sys
import os
import traceback
import csv
import pandas
import numpy
import matplotlib.pyplot as plt
from fisher_py import RawFile
from fisher_py.data.business import TraceType
from scipy.ndimage import gaussian_filter1d
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

    smoothed = gaussian_filter1d(i, sigma=5)

    plt.figure()
    graph = plt.plot(rt, i, label='Original', color='lightgray')
    graph2 = plt.plot(rt, smoothed, label='Smoothed', color='blue')

    # Gets the peak so that we know which one to integrate
    peak_index = numpy.argmax(smoothed)
    peak_rt= smoothed[peak_index]

    deriv = numpy.gradient(smoothed)

    # Calculates the left boundary by finding where you get back to baseline before the peak
    intensity = numpy.argmax(i)
    leftbound = peak_index
    while (intensity > 0 and intensity-smoothed[leftbound+1] < 0):
        leftbound -= 1
        intensity = smoothed[leftbound]

    # Calculates the right boundary by finding where you get back to baseline before the peak
    intensity = numpy.argmax(smoothed)
    rightbound = peak_index
    while (intensity > 0 and intensity-smoothed[rightbound-1] < 0):
        rightbound += 1
        intensity = smoothed[rightbound]

    # Makes a slice of the retention time and intensity lists so that it only involves the peak's curve
    simulated_baseline = numpy.linspace(smoothed[leftbound], smoothed[rightbound], num=(rightbound - leftbound))

    peakslicex = rt[leftbound:rightbound]
    peakslicey = i[leftbound:rightbound]
    peakslicey -= simulated_baseline
    peakslicey = numpy.clip(peakslicey, a_min=0, a_max=None)

    # Calculates the MA value (area under the peak's curve)
    MAval = trapezoid(peakslicey, peakslicex) * 60

    plt.scatter(rt[leftbound], i[leftbound])
    plt.scatter(rt[rightbound], i[rightbound])
    plt.savefig('myplot1.png')
    MAlist.append(MAval)
    print(f"Lysate MA = {MAval}, Lysate % acc = ", (MAval-2386741)/2386741*100)
    print(f"Media MA = {MAval}, Media % acc = ", (MAval-587693485)/587693485*100)

def filepather():
    dirinput = input("What's the name of the target directory? ")
    outputwrite = input("What file do you want to output? ")

    #path = Path.cwd()
    q = Path(dirinput)

    pathlist = []
    samplenames = []

    for child in q.iterdir():
        if child.name.endswith(".raw"):
            pathlist.append(child)

    for p in pathlist:
        samplenames.append(p.stem)
        extractor(str(p))

        # Saves the MAlist to a csv file that you specify
    with open(outputwrite, 'a', newline='') as csvfile:
        dawriter = csv.writer(csvfile)
        dawriter.writerow(samplenames)
        dawriter.writerow(MAlist)

filepather()
