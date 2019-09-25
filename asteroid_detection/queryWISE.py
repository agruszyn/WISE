import urllib.request as urll
from astropy.io import fits
import astropy.table
import numpy as np
import matplotlib.pyplot as plt
import requests
import time
import atpy
#t = atpy(table)
# # CALTECH code to build URL
# Examining the 'description' and 'catname' columns of this file,
# we find that the AllWISE Source Catalog has a cat name of "allwise_p3as_psd".
# https://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html
# can search by radial ascention
params = {'coadd_id': '0293p696_ac51',
          'band': 1,
          }
params['coaddgrp'] = params['coadd_id'][:2]
params['coadd_ra'] = params['coadd_id'][:4]
path = str.format(
    '{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits',
    **params)

#search for RA/Dec as functin of time

rimage = '03/0390/0390p605_ac51/0390p605_ac51-w1-int-3.fits?center=39.000655,60.577778&size=200pix'
image = '?center=39.000655,60.577778&size=200pix'
url = 'https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + rimage
print('https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + path)
#Sample Image Query:
#/ibe/data/wise/allsky/4band_p1bm_frm/6a/02206a/149/02206a149-w1-unc-1b.fits.gz?center=300,300pix&size=2arcmin&gzip=false
# /ibe/data/wise/allwise/p3am_cdd/03/0390/0390p605_ac51/0390p605_ac51-w1-int-3.fits?center=39.000655,60.577778&size=200pix

response = urll.urlopen('https://irsa.ipac.caltech.edu/ibe/data/wise/allsky/4band_p1bm_frm/6a/02206a/149/02206a149-w1-int-1b.fits?center=70,20&size=200pix')
html = response.read()
#tps://irsa.ipac.caltech.edu/ibe/data/wise/allsky/4band_p1bm_frm/6a/02206a/')
with open('images/test.fits', 'wb') as f:
    f.write(html)

fits_file = 'images/test.fits'
hdul = fits.open(fits_file)
data = hdul[0].data
plt.imshow(np.log10(data))
plt.show()

catalog = requests.get("https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mode=ascii")
with open('catalog.tbl', 'wb') as file:
    file.write(catalog.content)


def findfiles(html):
    record = False
    files = []
    word = []
    location = []
    string = ''
    time.sleep(0.25)
    r = requests.get(html)
    for i in range(len(r.text)):
        if r.text[(i - 6):i] == 'href="' and r.text[i] != '?' and r.text[i] != '/':
            record = True
        elif record and (r.text[i] == '"' or r.text[i] == "/"):
            record = False
            word = (string.join(word))
            location.append(word)
            word = []
        if record:
            word.append(r.text[i])

    for a in location:
        if '.' in a and '.gz' not in a and '.md5' not in a:
            files.append(a)
            print(files)
        elif '.' not in a and '106' in a:
            files.append(findfiles(html + '/' + a))
            print(files)
    return files


fits = []
tbl = []
dir = findfiles('https://irsa.ipac.caltech.edu/ibe/data/wise/allsky/4band_p1bm_frm/6a/02206a')
print('dir: ', dir[0])
for i in dir[0]:
    if '.fits' in i:
        fits.append(i)
    elif '.tbl' in i:
        tbl.append(i)

print('fits files: ', fits)
print('tbl files: ', tbl)

with open('fits_files.csv', 'w') as csvFile:
    for name in fits:
        csvFile.write(name)
        csvFile.write('\n')
csvFile.close()

with open('tbl_files.csv', 'w') as csvFile:
    for name in tbl:
        csvFile.write(name)
        csvFile.write('\n')
csvFile.close()

