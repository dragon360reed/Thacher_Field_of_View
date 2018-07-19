import healpy as hp
import numpy as np
import json
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
import csv
from scipy.special import erf
import astroquery as ast
from tqdm import tqdm
import pandas as pd

#Make pandas dataframe
csvinput = pd.read_csv("Thacher_Galaxies_sorted_by_dist.csv")
csvinput['used'] = np.chararray(len(csvinput))
csvinput['used'][:] = False
csvinput['Skycoord'] = np.nan
for galaxy in tqdm(range(len(csvinput))):
    csvinput['Skycoord'][galaxy] = coord.SkyCoord(csvinput['RA'][galaxy],csvinput['Dec'][galaxy],unit=(u.hourangle, u.deg),frame="fk5")

#The following outputs the distance between each galaxy and all other galaxies, takes a few minutes to run, Definitely not the most efficient way to do it
"""
Dis = {} #huge dictionary of all the distances
for galaxy in tqdm(range(len(csvinput))):
    Dis[csvinput['Name'][galaxy]] = {}
    for galaxy_compare in range(len(csvinput)):
        sep = csvinput['Skycoord'][galaxy].separation(csvinput['Skycoord'][galaxy_compare])
        Dis[csvinput['Name'][galaxy]][csvinput['Name'][galaxy_compare]] = sep
"""
#Read results from block above from pickle file
Dis = pickle.load(open('Distances_dict.pkl','rb'))

#Make dictionary to np array, so it's easier to manipulate
Dis_array = np.array([[Dis[j][i].arcminute for i in csvinput['Name']] for j in csvinput['Name']])
for i in tqdm(range(len(csvinput))):
    for j in range(len(csvinput)):
        Darray[i][j]=Dis[csvinput['Name'][i]][csvinput['Name'][j]].deg

#Find galaxy pairs that are within our field of view: 20.7725 x 20.7456 arcmin
within_field = []
for i in range(430):
    for j in range(430):
        if Darray[i][j] <= 0.33 and Darray[i][j]!=0.:
            within_field.append([i, j])

#Plot the galaxies and the red lines indicate the ones that fit within one frame:
plt.ion()
plt.figure(figsize=(10,5))
plt.scatter([i.ra.deg for i in csvinput['Skycoord'].values],[i.dec.deg for i in csvinput['Skycoord'].values])
plt.title('Thacher Galaxies Catalogue')
plt.xlabel('RA')
plt.ylabel('Dec')
for i in within_field:
    plt.plot([csvinput['Skycoord'][i[0]].ra.deg,csvinput['Skycoord'][i[1]].ra.deg],[csvinput['Skycoord'][i[0]].dec.deg,csvinput['Skycoord'][i[1]].dec.deg],'r-')


"""
Some modification of the translation of C++ code, not sure if it is exactly what Draco wants

csvinput['Within_field'] = np.chararray(len(csvinput))
csvinput['Within_field'][:] = False

r = 18. #I need to check on this, in arcminutes

one_field = [] #Forms an empty list to put galaxies within one field of view in

for row in range(len(csvinput)):

    if csvinput['used'][row] == False:

            ra1 = np.cos(csvinput['Dec'][row])*csvinput['RA'][row] #RA requires a cosine factor.
            dec1 = csvinput['Dec'][row]

            ra2 = np.cos(csvinput['Dec'][row+1])*csvinput['RA'][row+1]#RA requires a cosine factor.                dec2 = (row+1)[2]

            dist = np.sqrt(dist = np.sqrt((ra2-ra1)**2. + (dec2-dec1)**2.))

            if dist < r:
                    one_field.append(dist)
                    csvinput['Within_field'][row] = True
            else:
                print ("Cannot fit in frame")

            csvinput['used'][row] == True

    elif csvinput['used'][row] == True:

            row = row+1 #We need to skip the ones which have already been chosen

csvinput.to_csv("Thacher_Galaxies_within_field_query.csv")
"""
