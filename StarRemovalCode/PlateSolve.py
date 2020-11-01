from StarRemovalCode.Client import Client
import os
import time
import astropy as ap
from astropy.io import fits

from urllib.parse import urlparse, urlencode, quote
from urllib.request import urlopen, Request
from urllib.error import HTTPError

def getStars(fileIn, showProgress):
    server = Client.default_url

    c = Client(Client.default_url)
    c.login("ujazexnqbmqqikyl")

    kwargs = dict()
    #kwargs["allow_commercial_use"]=False
    #kwargs["allow_modifications"]=False
    #kwargs["publicly_visible"]=False
    kwargs['scale_type'] = "ul"
    #kwargs['scale_lower'] = 0.1
    #kwargs['scale_upper'] = 60.0
    #kwargs['scale_units'] = "arcsecperpix"
    #kwargs['downsample_factor'] = 2
    #kwargs['tweak_order'] = 2
    #kwargs['use_sextractor'] = False
    #kwargs['positional_error'] = 1
    #kwargs['crpix_center'] = True



    upres = c.upload(fileIn, **kwargs)

    if showProgress:
        print("Uploading...")

    stat = upres['status']
    if stat != 'success':
        print('Upload failed: status', stat)
        print(upres)
        sys.exit(-1)
    else:
        if showProgress:
            print("Upload successful\n")
            print("Your image is in the Astrometry.net queue. This may take several minutes.\n")

    sub_id = upres['subid']
    solved_id = None
    success = False

    while True:
        stat = c.sub_status(sub_id, justdict=True)
        jobs = stat.get('jobs', [])
        if len(jobs):
            for j in jobs:
                if j is not None:
                    break
            if j is not None:
                solved_id = j
                break
        time.sleep(1)

    if showProgress:
        print("Platesolving...")

    while True:
        stat = c.job_status(solved_id, justdict=True)
        if stat.get('status','') in ['success']:
            success = (stat['status'] == 'success')
            break
        time.sleep(1)

    if showProgress:
        print("Image successfully platesolved\n")
    if success:
        url = server.replace('/api/', '/wcs_file/%i/' % solved_id)#server.replace('/api/', '/corr_file/%i/' % solved_id)
        f = urlopen(url)
        pos = fits.open(f)
        
        return pos[0].header

