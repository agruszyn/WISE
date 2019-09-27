import urllib.request as urll
from astropy.io import fits
import astropy.table
import os
import numpy as np
import matplotlib.pyplot as plt
import requests
import time
import atpy


params = {'scan_id': '50077b',
          'frame_num': 150,
          'band': 1,
          }

params['scangrp'] = params['scan_id'][-2:]
path = str.format(
    '{scangrp:s}/{scan_id:s}/{frame_num:03d}/{scan_id:s}{frame_num:03d}-w{band:1d}-int-1b.fits',
    **params)
neowiser = 'https://irsa.ipac.caltech.edu/ibe/data/wise/neowiser/p1bm_frm/' + path

radecQuery = 'https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?POS=0,0'#'https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?<parameters>'

URL = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?'
spatialQuery = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=box&catalog=allwise_p3as_psd&objstr=00h+42m+44.32s+41d+16m+08.5s&size=300&outfmt=1'
catalog = 'name'
spatial = 'cone' #box, polygon, upload, none
radius = '0.5'
objstr = 3# 00h+42m+42.32s+41d+16m+08.5s
radunits = 'arcsec' # arcmin, deg
outfmt = 0 #0:HTML, 1: ASCII, 2: SVC, 3: VO Table, 6: XML

#html = requests.get("https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?outfmt=1&searchForm=MO&spatial=cone&catalog=allsky_4band_p1bs_psd&moradius=5&mobj=smo&mobjstr=324")
html = requests.get(neowiser)
file = open('images/R0D0image.fits', 'wb')
file.write(html.content)

fits_file = 'images/R0D0image.fits'
hdul = fits.open(fits_file)
data = hdul[0].data
plt.imshow(np.log10(data))
plt.show()
file.close()
#print(os.getcwd())
#tab = atpy.Table('R0D0metadata.tbl')
#print(tab.columns)
