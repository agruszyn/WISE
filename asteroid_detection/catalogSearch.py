import urllib.request as urll
from astropy.io import fits
import astropy.table
import os
import numpy as np
import matplotlib.pyplot as plt
import requests
import time
import atpy

URL = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?'
URL2 = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=box&catalog=allwise_p3as_psd&objstr=00h+42m+44.32s+41d+16m+08.5s&size=300&outfmt=1'
catalog = 'name'
spatial = 'cone' #box, polygon, upload, none
radius = '0.5'
objstr = 3# 00h+42m+42.32s+41d+16m+08.5s
radunits = 'arcsec' # arcmin, deg
outfmt = 0 #0:HTML, 1: ASCII, 2: SVC, 3: VO Table, 6: XML

#html = requests.get("https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?outfmt=1&searchForm=MO&spatial=cone&catalog=allsky_4band_p1bs_psd&moradius=5&mobj=smo&mobjstr=324")
html = requests.get(URL2)
file = open('coneSearch.tbl', 'wb')
file.write(html.content)
file.close()
print(os.getcwd())
tab = atpy.Table('coneSearch.tbl')
print(tab.columns)
