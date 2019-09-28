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


def get_asteroid_sources(catname, asteroid):
    if not os.path.isfile(asteroid + '_Catalog.tbl'):
        html = requests.get(URL + 'cgi-bin/Gator/nph-query?outfmt=1&searchForm=MO&spatial=cone&catalog=' + catname + '&moradius=5&mobj=smo&mobjstr=' + asteroid)
        file = open(asteroid + '_Catalog.tbl', 'wb')
        file.write(html.content)
        file.close()
    return atpy.Table(asteroid + '_Catalog.tbl')


def pick_sourcepoint(catalog):
    # get size
    length = np.size(catalog.mjd)
    # pick sourcepoint
    return np.random.randint(low=0, high=length)


def get_datetime(catalog, index):
    t = Time(catalog.mjd[index], format='mjd')
    t.format = 'fits'
    t = str(t).replace('T', ' ', 1)
    return t


def get_positional_images(catalog, index):
    search = 'https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?POS=' + str(catalog.ra[index]) + ',' + str(catalog.dec[index])
    html = requests.get(search)
    file = open('fitsData.tbl', 'wb')
    file.write(html.content)
    file.close()
    return atpy.Table('fitsData.tbl')


def find_image(catalog, date):
    for i in range(len(catalog.date_obs)):
        if utc in catalog.date_obs[i]:
            return {'scan_id': str(catalog.scan_id[i]),
                    'frame_num': catalog.frame_num[i],
                    'band': catalog.band[i],
                    }


def download_image(asteroid, params):
    image = 'images/' + asteroid + '_' + str.format('{scan_id:s}{frame_num:03d}-w{band:1d}', **params) + '-int-1b.fits'

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
        return image


def display_fits(fits):
    fits_file = fits
    hdul = fits.open(fits_file)
    data = hdul[0].data
    hdul.close()
    plt.imshow(np.log10(data))
    plt.show()
    return data


def filter(data):
    for a in range(len(data)):
        for b in range(len(data[a])):
            if np.isnan(data[a, b]):
                data[a, b] = 255
    return data

#OBJECT = 'Ceres'
OBJECT = 'Gletorrence'
URL = 'https://irsa.ipac.caltech.edu/'
catnames = ['neowiser_p1bs_psd', 'neowiser_p1ba_mch', 'neowiser_p1bs_frm', 'neowiser_p1bl_lod']
catalogName = 'neowiser_p1bs_psd'

asteroidMetadata = get_asteroid_sources(catalogName, OBJECT)

element = pick_sourcepoint(asteroidMetadata)

utc = get_datetime(asteroidMetadata, element)

coordinateMetadata = get_positional_images(asteroidMetadata, element)

scan_id, frame_num, band = find_image(coordinateMetadata, utc)

params = find_image(coordinateMetadata, utc)

image = download_image(OBJECT, params)

w = wcs.WCS(image)
x, y = w.wcs_world2pix(asteroidMetadata.ra[element], asteroidMetadata.dec[element], 0)

X = display_fits(image)

X = filter(X)



fx = dip.fftshift(dip.fft2(X))
radius = 30
center = len(fx)/2

f = np.log(np.abs(fx))
for a in range(len(f)):
    for b in range(len(f[a])):
        if np.square(center - a) + np.square(center - b) < np.square(radius):
            f[a, b] = 0
plt.imshow(np.log10(x))
plt.show()

plt.imshow(np.log(np.abs(dip.ifft2(fx))))
plt.show()