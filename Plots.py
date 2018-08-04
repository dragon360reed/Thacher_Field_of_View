from matplotlib.patches import Rectangle, Polygon
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import astropy.units as u
import astropy.coordinates as coord


Final_df_avcut_hstcut = pd.read_csv('Final_df_avcut_hstcut.csv')
#Dis_array = pickle.load(open("Dis_array.pkl","rb"))
within_field = pickle.load(open("within_field.pkl","rb"))

fov = 0.25


SkyCoord = coord.SkyCoord(Final_df_avcut_hstcut['RA'].values,Final_df_avcut_hstcut['Dec'].values,unit=(u.deg, u.deg))
ra_values = SkyCoord.ra.degree
dec_values = SkyCoord.dec.degree

Pointings = {}
Pointings['RA'] = []
Pointings['Dec'] = []
Pointings["Galaxies"] = []

#Singles:
singles = []
for g in [i for i in range(len(Final_df_avcut_hstcut)) if i not in np.unique(within_field)]:
    singles.append(g)
    ra = Final_df_avcut_hstcut['RA'][g]
    dec = Final_df_avcut_hstcut['Dec'][g]
    Pointings['RA'].append(ra)
    Pointings['Dec'].append(dec)
    Pointings['Galaxies'].append([Final_df_avcut_hstcut['GWGC_name'][g]])
    #rect = Rectangle((ra-fov/2,dec-fov/2),fov,fov,linewidth=1,edgecolor='orange',facecolor='none')
    #ax.add_patch(rect)

def centroid(array):
    length = array.shape[0]
    sum_x = np.sum(array[:, 0])
    sum_y = np.sum(array[:, 1])
    return sum_x/length, sum_y/length

#Triangles
triangles = []
for i in range(len(within_field)-1):
    if within_field[i][0] == within_field[i+1][0] and [within_field[i][1],within_field[i+1][1]] in within_field:
        a = within_field[i][0]
        b = within_field[i][1]
        c = within_field[i+1][1]
        if [list(np.array(within_field).T[0]).count(a) == 2, list(np.array(within_field).T[0]).count(b) == 2, list(np.array(within_field).T[0]).count(c) ==2].count(True) >=2:
            triangles.append([a,b,c])

#Pairs:
pairs = []
for i in within_field:
    if list(np.array(within_field).T[0]).count(i[0]) ==1 or list(np.array(within_field).T[0]).count(i[1]) ==1 :
        pairs.append(i)

for t in triangles:
    array = np.array([[Final_df_avcut_hstcut['RA'][t[0]],Final_df_avcut_hstcut['Dec'][t[0]]],
    [Final_df_avcut_hstcut['RA'][t[1]],Final_df_avcut_hstcut['Dec'][t[1]]],
    [Final_df_avcut_hstcut['RA'][t[2]],Final_df_avcut_hstcut['Dec'][t[2]]]])
    cra, cdec = centroid(array)
    #rect = Rectangle((cra-fov/2,cdec-fov/2),fov,fov,linewidth=1,edgecolor='purple',facecolor='none')
    Pointings['RA'].append(cra)
    Pointings['Dec'].append(cdec)
    Pointings['Galaxies'].append([Final_df_avcut_hstcut['GWGC_name'][t[0]],Final_df_avcut_hstcut['GWGC_name'][t[1]],Final_df_avcut_hstcut['GWGC_name'][t[2]]])
    #ax.add_patch(rect)

for p in pairs:
    cra = (Final_df_avcut_hstcut['RA'][p[0]] + Final_df_avcut_hstcut['RA'][p[1]])/2.
    cdec = (Final_df_avcut_hstcut['Dec'][p[0]] + Final_df_avcut_hstcut['Dec'][p[1]])/2.
    #rect = Rectangle((cra-fov/2,cdec-fov/2),fov,fov,linewidth=1,edgecolor='green',facecolor='none')
    Pointings['RA'].append(cra)
    Pointings['Dec'].append(cdec)
    Pointings['Galaxies'].append([Final_df_avcut_hstcut['GWGC_name'][p[0]],Final_df_avcut_hstcut['GWGC_name'][p[1]]])
    #ax.add_patch(rect)

#Clusters:
clusters = [i for i in range(len(Final_df_avcut_hstcut)) if i not in np.unique(singles) and i not in np.unique(triangles) and i not in np.unique(pairs)]
for i in clusters:
    plt.annotate(str(i),(Final_df_avcut_hstcut['RA'][i],Final_df_avcut_hstcut['Dec'][i]))

manual = [[190,208,418,978],[978,259],[923,530,205],[236,552],[542,924,316],[332,659,924,316]]
for m in manual:
    array = np.array([[Final_df_avcut_hstcut['RA'][i],Final_df_avcut_hstcut['Dec'][i]] for i in m])
    cra, cdec = centroid(array)
    #rect = Rectangle((cra-fov/2,cdec-fov/2),fov,fov,linewidth=1,edgecolor='black',facecolor='none')
    Pointings['RA'].append(cra)
    Pointings['Dec'].append(cdec)
    Pointings['Galaxies'].append([Final_df_avcut_hstcut['GWGC_name'][i] for i in m])
    #ax.add_patch(rect)

#for g in np.unique(within_field):
#    plt.annotate(Final_df_avcut_hstcut['GWGC_name'][g],(Final_df_avcut_hstcut['RA'][g],Final_df_avcut_hstcut['Dec'][g]))
Pointings = pd.DataFrame.from_dict(Pointings)
Pointings.to_csv("Pointings.csv")

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(ra_values,dec_values)
plt.title('Thacher Galaxies Catalogue')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.subplots_adjust(top=0.95,bottom=0.0)
for i in within_field:
    plt.plot([ra_values[i[0]], ra_values[i[1]]],[dec_values[i[0]],dec_values[i[1]]],'r-')
Pointings_coord = coord.SkyCoord(Pointings['RA'],Pointings['Dec'],unit=(u.deg, u.deg))
pointings_ra = Pointings_coord.ra.degree
pointings_dec = Pointings_coord.dec.degree
for p in range(len(Pointings['RA'])):
    r = pointings_ra[p]
    d = pointings_dec[p]
    ra_offset = fov/2/np.cos(np.radians(d))
    dec_offset = fov/2
    x = [r-ra_offset,r+ra_offset,r+ra_offset,r-ra_offset]
    y = [d+dec_offset,d+dec_offset,d-dec_offset,d-dec_offset]
    ax.add_patch(Polygon(xy=list(zip(x,y)), fill=False, edgecolor='orange'))
plt.grid(True)
