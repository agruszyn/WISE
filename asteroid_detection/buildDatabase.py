import urllib.request as urll
from astropy.io import fits
# import astropy.table
# import astropy.time
import io
from astropy.time import Time
import time
import dippykit as dip
import os
import numpy as np
import matplotlib.pyplot as plt
import requests
import atpy
import json
from astropy import wcs
#import pyvo


def get_asteroid_sources(catname, asteroid, file):
    print('get sources')
    html = requests.get(URL + 'cgi-bin/Gator/nph-query?outfmt=1&searchForm=MO&spatial=cone&catalog=' + catname + '&moradius=5&mobj=smo&mobjstr=' + asteroid)
    file.write(html.content)
    file.close()
    return atpy.Table('sources.tbl')


def get_datetime(catalog, index):
    print('get datetime')
    t = Time(catalog.mjd[index], format='mjd')
    t.format = 'fits'
    t = str(t).replace('T', ' ', 1)
    return t


def find_image(catalog, date):
    print('find image')
    print(date)
    for i in range(len(catalog.date_obs)):
        if date in catalog.date_obs[i]:
            return {'scan_id': str(catalog.scan_id[i]),
                    'frame_num': str(catalog.frame_num[i]),
                    'band': str(catalog.band[i]),
                    }


def get_positional_images(catalog, index, utc, file):
    print('find images at position')
    search = 'https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?POS=' + str(catalog.ra[index]) + ',' + str(catalog.dec[index])
    html = requests.get(search)
    file.write(html.content)
    file.close()
    return find_image(atpy.Table('positionTable.tbl'), utc)



OBJECT = 'Ceres'
URL = 'https://irsa.ipac.caltech.edu/'
catalogName = 'neowiser_p1bs_psd'

filecontent = open('sources.tbl', 'wb')
asteroidMetadata = get_asteroid_sources(catalogName, OBJECT, filecontent)

ceres = []
dictionary = {}

#for i in range(len(asteroidMetadata)):
for i in range(5):
    newfile = open('positionTable.tbl', 'wb')
    utc = (get_datetime(asteroidMetadata, i))
    params = get_positional_images(asteroidMetadata, i, utc, newfile)
    if params is not None:
        key = params['scan_id'] + params['frame_num'] + params['band']
        if key not in dictionary:
            dictionary[key] = {}
        dictionary[key]['params'] = params
        if 'asteroids' not in dictionary[key]:
            dictionary[key]['asteroids'] = []
        asteroid = {}
        asteroid['ra'] = str(asteroidMetadata.ra[i])
        asteroid['dec'] = str(asteroidMetadata.dec[i])
        asteroid['name'] = OBJECT
        asteroid['date'] = utc
        dictionary[key]['asteroids'].append(asteroid)
        print('scans left: ', len(asteroidMetadata) - i)
        time.sleep(5)

file = open('mylist.json', 'w')
json.dump(dictionary, file)
file.close()
