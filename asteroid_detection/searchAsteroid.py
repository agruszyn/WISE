import urllib.request as urll
from astropy.io import fits
import astropy.table
import astropy.time
import os
import numpy as np
import matplotlib.pyplot as plt
import requests
import time
import atpy
import pyvo

if not os.path.isfile('neo0.tbl'):
    OBJECT = 'ceres'
    URL = 'https://irsa.ipac.caltech.edu/'
    catnames = ['neowiser_p1bs_psd', 'neowiser_p1ba_mch', 'neowiser_p1bs_frm', 'neowiser_p1bl_lod']
    #catname[0] is moving object enables, the others are not
    html = requests.get(URL + 'cgi-bin/Gator/nph-query?outfmt=1&searchForm=MO&spatial=cone&catalog=' + catnames[0] + '&moradius=5&mobj=smo&mobjstr=324')
    file = open('neo0.tbl', 'wb')
    file.write(html.content)
    file.close()


tab = atpy.Table('neo0.tbl')
#print(tab.dec)
#for column name info:
#http://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec2_1a.html
#http://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec2_1a.html
#/ibe/search/wise/neowiser/p1bm_frm?POS=0,0

t = astropy.time.Time(tab.mjd[2], format='mjd')
t.format = 'fits'
print('HERE', t)

search = 'https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?POS=' + str(tab.ra[2]) + ',' + str(tab.dec[2])



html2 = requests.get(search)
file = open('fitsData.tbl', 'wb')
file.write(html2.content)
file.close()

fitsData = atpy.Table('fitsData.tbl')
print(fitsData.date_obs)
utc = str(t).replace('T', ' ', 1)
print('utc: ', utc)


result = np.where('2013-12-25' in fitsData.date_obs)
for i in range(len(fitsData.date_obs)):
    if utc in fitsData.date_obs[i]:
        print(fitsData.scan_id[i], ' ', fitsData.frame_num[i], ' ', fitsData.band[i])

params = {'scan_id': '44558a',
          'frame_num': 174,
          'band': 1,
          }
if not os.path.isfile('images/ceres.tbl'):
    params['scangrp'] = params['scan_id'][-2:]
    path = str.format(
        '{scangrp:s}/{scan_id:s}/{frame_num:03d}/{scan_id:s}{frame_num:03d}-w{band:1d}-int-1b.fits',
        **params)
    url = 'https://irsa.ipac.caltech.edu/ibe/data/wise/neowiser/p1bm_frm/' + path

    r = requests.get(url)
    file = open('images/ceres.fits', 'wb')
    file.write(r.content)
    file.close()

fits_file = 'images/ceres.fits'
hdul = fits.open(fits_file)
data = hdul[0].data
plt.imshow(np.log10(data))
plt.show()

