import urllib.request as urll
from astropy.io import fits
# import astropy.table
# import astropy.time
from astropy.time import Time
import dippykit as dip
import os
import numpy as np
import matplotlib.pyplot as plt
import requests
import atpy
from astropy import wcs
#import pyvo

#OBJECT = 'ceres'
OBJECT = 'gletorrence'
URL = 'https://irsa.ipac.caltech.edu/'
catnames = ['neowiser_p1bs_psd', 'neowiser_p1ba_mch', 'neowiser_p1bs_frm', 'neowiser_p1bl_lod']

if not os.path.isfile('neo0.tbl'):
    html = requests.get(URL + 'cgi-bin/Gator/nph-query?outfmt=1&searchForm=MO&spatial=cone&catalog=' + catnames[0] + '&moradius=5&mobj=smo&mobjstr=324')
    file = open('neo0.tbl', 'wb')
    file.write(html.content)
    file.close()



tab = atpy.Table('neo0.tbl')
count = np.size(tab.mjd)
index = np.random.randint(low=0, high=count)
print('index: ', index)

t = Time(tab.mjd[index], format='mjd')
t.format = 'fits'

search = 'https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?POS=' + str(tab.ra[index]) + ',' + str(tab.dec[index])
print(tab.ra[index], ' and ', tab.dec[index])

#if not os.path.isfile('fitsData.tbl'):
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
print('UTC: ', utc)
print(fitsData.date_obs)
params = {'scan_id': scan_id[0],
          'frame_num': frame_num[0],
          'band': band[-1],
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

w = wcs.WCS(image)
x, y = w.wcs_world2pix(tab.ra[index], tab.dec[index], 0)
print('lon: ', x, 'lat   : ', y)

fits_file = image
hdul = fits.open(fits_file)
data = hdul[0].data
header = hdul[0].header
hdul.close()
print(type(data))
plt.imshow(np.log10(data))
plt.show()

for a in range(len(data)):
    for b in range(len(data[a])):
        if np.isnan(data[a,b]):
            data[a, b] = 255
            print(a, b)

fx = dip.fftshift(dip.fft2(data))
radius = 30
center = len(fx)/2

f = np.log(np.abs(fx))
for a in range(len(f)):
    for b in range(len(f[a])):
        if np.square(center - a) + np.square(center - b) < np.square(radius):
            f[a, b] = 0
plt.imshow(np.log10(data))
plt.show()

plt.imshow(np.log(np.abs(dip.ifft2(fx))))
plt.show()

# CERES
#342 x
#925 y
# CRPIX1  =                508.5 / Reference X pixel
# CRPIX2  =                508.5 / Reference Y pixel
# CRVAL1  =    0.433682830922963 / [deg] Image center RA
# CRVAL2  =     14.0325824797345 / [deg] Image center Dec
# CTYPE1  = 'RA---SIN-SIP'       / Sin projection with SIP coefficients
# CTYPE2  = 'DEC--SIN-SIP'       / Sin projection with SIP coefficients
# CD1_1   = 0.000698823672152697 / WCS rotation matrix element
# CD1_2   = 0.000312769537857026 / WCS rotation matrix element
# CD2_1   = 0.000314830254725055 / WCS rotation matrix element
# CD2_2   = -0.000694249531937869 / WCS rotation matrix element
# WCROTA2 =     204.252259008431 / [deg] CCW rotation of RA at CRPIX1,2
# PA      =     155.747740991569 / [deg] Rotation of +Y EofN at CRPIX1,2
# WCDELT1 = -0.000766467621006409 / [deg/pix] X-axis scale
# WCDELT2 = 0.000761450718305031 / [deg/pix] Y-axis scale
# PXSCAL1 =    -2.75928343562307 / [arcsec/pixel] X-axis scale
# PXSCAL2 =     2.74122258589811 / [arcsec/pixel] Y-axis scale