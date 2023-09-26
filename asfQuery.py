"""
Created on Tue Jun 14 11:26:10 2022

asfQuery.py

path: str '64'
frame: best to leave this as None unless there's a specific frame you exclusively want
start: str '2010-05-01T00:00:00Z'
end: str '2010-05-01T00:00:00Z'

You can give a wkt polygon, bounds, or a point. It will prefer them in that order.
polygon: Must be in WKT format already: 'POLYGON((-160.1985 21.7373,-159.2377 21.8944,-159.2195 22.1633,-159.3706 22.2807,-159.7513 22.2192,-160.2711 21.9392,-160.2952 21.8046,-160.1985 21.7373))'
bounds: 'S, N, W, E'
point: 'lon, lat'

Ideally, you should put the polygon and the bounds in localParams.  This script
will use the polygon to search for only the intersecting frames, and then isce
stack processor will use the bounds to only process the bursts intersecting
those bounds.

@author: km
"""
import requests, pandas as pd
from lxml import html
import shapely, shapely.geometry, shapely.wkt
import sys
import time

def getGran(path, start, end, sat, bounds, poly):

    # path=ps.path
    # start=ps.start
    # end=ps.end
    # sat=ps.sat
    # bounds=ps.bounds
    # point=ps.point
    # poly=ps.poly

    fmt = 'CSV'
    flightDirection = None
    baseurl = 'https://api.daac.asf.alaska.edu/services/search/param'
    data = dict(platform=sat,
      processingLevel='SLC',
      output=fmt)


    if poly:
        print('Using polygon')
        data['intersectsWith'] = poly
    elif bounds:
        print('Using bounds.  Usually it is better to use a polygon instead.')
        miny, maxy, minx, maxx = bounds.split(sep=',')
        roi = shapely.geometry.box(float(minx), float(miny), float(maxx), float(maxy))
        polygonWKT = roi.wkt
        data['intersectsWith'] = polygonWKT
    else:
        print('Need to specify the polygon search area')
        sys.exit(1)

    if path:
        data['relativeOrbit'] = path
    if start:
        data['start'] = start
    if end:
        data['end'] = end
    if flightDirection:
        data['flightDirection'] = flightDirection

    print('Searching for data')
    r = requests.get(baseurl, params=data, timeout=100)
    print(r)
    with open('out.csv', 'w') as (j):
        j.write(r.text)
    slcUrls = pd.read_csv('out.csv')['URL']
    sizesMB = pd.read_csv('out.csv')['Size (MB)']
    dates = pd.read_csv('out.csv')["Acquisition Date"]
    print('Found ' + str(round((sum(sizesMB)/1000),2)) + ' Gb of data')
    gran = pd.read_csv('out.csv')['Granule Name']
    
    return (slcUrls, gran,dates,r)


def get_with_retry(url, max_retries=3):
    retries = 0
    while retries < max_retries:
        response = requests.get(url, stream=True)
        if response.status_code != 429:
            return response
        wait_time = (2 ** retries)
        print(f"Rate limited! Waiting for {wait_time} seconds before retrying...")
        time.sleep(wait_time)
        retries += 1
    response.raise_for_status()


def get_orbit_url(granuleName):
    """Retrieve precise orbit file for a specific Sentinel-1 granule.
    Precise orbits available ~3 weeks after aquisition.
    Parameters
    ----------
    granuleName : str
        ASF granule name, e.g.:
        S1B_IW_SLC__1SDV_20171117T015310_20171117T015337_008315_00EB6C_40CA
    url : str
        website with simple list of orbit file links
    Returns
    -------
    orbitUrl :  str
        url pointing to matched orbit file
    """
    
    urlPrecise='https://s1qc.asf.alaska.edu/aux_poeorb'
    urlResorb='https://s1qc.asf.alaska.edu/aux_resorb'

    sat = granuleName[:3]
    date = granuleName[3:]
    print(f"retrieving precise orbit URL for {sat}, {date}")
    
    try:
        r = get_with_retry(urlPrecise)
        webpage = html.fromstring(r.content)
        orbits = webpage.xpath('//a/@href')
        df = pd.DataFrame(dict(orbit=orbits))
        dfSat = df[df.orbit.str.startswith(sat)].copy()
        dayBefore = pd.to_datetime(date) - pd.to_timedelta(1, unit='d')
        dayBeforeStr = dayBefore.strftime('%Y%m%d')
        dfSat.loc[:, 'startTime'] = dfSat.orbit.str[42:50]
        match = dfSat.loc[(dfSat.startTime == dayBeforeStr, 'orbit')].values[0]
        orbitUrl = f"{urlPrecise}/{match}"
    except:
        try:
            print('using resorb for ' + granuleName +' (it is probably too recent)')
            r = get_with_retry(urlResorb)
            webpage = html.fromstring(r.content)
            orbits = webpage.xpath('//a/@href')
            df = pd.DataFrame(dict(orbit=orbits))
            dfSat = df[df.orbit.str.startswith(sat)].copy()
            dayBefore = pd.to_datetime(date) - pd.to_timedelta(1, unit='d')
            dayBeforeStr = dayBefore.strftime('%Y%m%d')
            dfSat.loc[:, 'startTime'] = dfSat.orbit.str[42:50]
            match = dfSat.loc[(dfSat.startTime == dayBeforeStr, 'orbit')].values[0]
            orbitUrl = f"{urlResorb}/{match}"
            # os.system('mv ./orbits/2*/*EOF ./orbits/')
        except Exception as e:
            print(f"Error encountered: {e}")
            raise RuntimeError("Both precise and resorb URL retrieval failed!") from e

    return orbitUrl
