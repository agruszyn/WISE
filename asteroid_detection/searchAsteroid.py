import urllib.request as urll
from astropy.io import fits
# import astropy.table
# import astropy.time
from astropy.time import Time
import os
import numpy as np
import matplotlib.pyplot as plt
import requests
import atpy
#import pyvo

OBJECT = 'ceres'
URL = 'https://irsa.ipac.caltech.edu/'
catnames = ['neowiser_p1bs_psd', 'neowiser_p1ba_mch', 'neowiser_p1bs_frm', 'neowiser_p1bl_lod']

if not os.path.isfile('neo0.tbl'):
    html = requests.get(URL + 'cgi-bin/Gator/nph-query?outfmt=1&searchForm=MO&spatial=cone&catalog=' + catnames[0] + '&moradius=5&mobj=smo&mobjstr=324')
    file = open('neo0.tbl', 'wb')
    file.write(html.content)
    file.close()


tab = atpy.Table('neo0.tbl')

t = Time(tab.mjd[2], format='mjd')
t.format = 'fits'

search = 'https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?POS=' + str(tab.ra[2]) + ',' + str(tab.dec[2])

if not os.path.isfile('fitsData.tbl'):
    html2 = requests.get(search)
    file = open('fitsData.tbl', 'wb')
    file.write(html2.content)
    file.close()

fitsData = atpy.Table('fitsData.tbl')
utc = str(t).replace('T', ' ', 1)

scan_id = []
frame_num = []
band = []
result = np.where(utc in fitsData.date_obs)

for i in range(len(fitsData.date_obs)):
    if utc in fitsData.date_obs[i]:
        #print(fitsData.scan_id[i], ' ', fitsData.frame_num[i], ' ', fitsData.band[i])
        scan_id.append(str(fitsData.scan_id[i]))
        frame_num.append(fitsData.frame_num[i])
        band.append(fitsData.band[i])

params = {'scan_id': scan_id[0],
          'frame_num': frame_num[0],
          'band': band[0],
          }
image = 'images/' + OBJECT + '_' + str.format('{scan_id:s}{frame_num:03d}-w{band:1d}', **params) + '-int-1b.fits'

if not os.path.isfile(image):
    params['scangrp'] = params['scan_id'][-2:]
    path = str.format(
        '{scangrp:s}/{scan_id:s}/{frame_num:03d}/{scan_id:s}{frame_num:03d}-w{band:1d}-int-1b.fits',
        **params)
    url = 'https://irsa.ipac.caltech.edu/ibe/data/wise/neowiser/p1bm_frm/' + path

    r = requests.get(url)
    file = open(image, 'wb')
    file.write(r.content)
    file.close()

fits_file = image
hdul = fits.open(fits_file)
data = hdul[0].data
header = hdul[0].header
print(header)
plt.imshow(np.log10(data))
plt.show()

