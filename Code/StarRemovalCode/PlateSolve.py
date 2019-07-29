import argparse
import Client
from Client import Client
import os
import time

from urllib.parse import urlparse, urlencode, quote
from urllib.request import urlopen, Request
from urllib.error import HTTPError
import json
python2json = json.dumps


ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required = True)
args = vars(ap.parse_args())

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



upres = c.upload(args["input"], **kwargs)

stat = upres['status']
if stat != 'success':
    print('Upload failed: status', stat)
    print(upres)
    sys.exit(-1)

sub_id = upres['subid']
solved_id = None
success = False

while True:
    stat = c.sub_status(sub_id, justdict=True)
    print('Got status:', stat)
    jobs = stat.get('jobs', [])
    if len(jobs):
        for j in jobs:
            if j is not None:
                break
        if j is not None:
            print('Selecting job id', j)
            solved_id = j
            break
    time.sleep(5)

while True:
    stat = c.job_status(solved_id, justdict=True)
    print('Got job status:', stat)
    if stat.get('status','') in ['success']:
        success = (stat['status'] == 'success')
        break
    time.sleep(5)


if success:
    # we have a jobId for retrieving results
    retrieveurls = []
    #url = server.replace('/api/', '/wcs_file/%i' % solved_id)

    #retrieveurls.append((url, True))
    #url = server.replace('/api/', '/kml_file/%i/' % solved_id)

    #retrieveurls.append((url, True))
    url = server.replace('/api/', '/new_fits_file/%i/' % solved_id)

    #retrieveurls.append((url, True))

    url = server.replace('/api/', '/corr_file/%i/' % solved_id)

    #retrieveurls.append((url, True))

    out = args["input"][0:-5]+"_corr"+args["input"][-5:]
    
    print('Retrieving file from', url, 'to', out)
    f = urlopen(url)
    txt = f.read()
    w = open(out, 'wb')
    w.write(txt)
    w.close()
    print('Wrote to', out)

    result = c.annotate_data(solved_id)
    print("annotations")
    print(result["annotations"])
    with open("results.json",'w') as f:
        f.write(python2json(result))

