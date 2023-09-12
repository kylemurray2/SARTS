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
import os, requests, pandas as pd
from lxml import html
import shapely, shapely.geometry, shapely.wkt
import sys
from datetime import timedelta

def getGran(path, start, end, sat, bounds, poly):

    # path=ps.path
    # start=ps.start
    # end=ps.end
    # sat=ps.sat
    # bounds=ps.bounds
    # point=ps.point
    # poly=ps.poly

    fmt = 'CSV'
    orbit = None
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
    date = granuleName[17:25]
    print(f"retrieving precise orbit URL for {sat}, {date}")
    
    try:
        r = requests.get(urlPrecise,stream=True)
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
        print('using resorb for this date (it is probably too recent)')
        r = requests.get(urlResorb)
        webpage = html.fromstring(r.content)
        orbits = webpage.xpath('//a/@href')
        df = pd.DataFrame(dict(orbit=orbits))
        dfSat = df[df.orbit.str.startswith(sat)].copy()
        dayBefore = pd.to_datetime(date) - pd.to_timedelta(1, unit='d')
        dayBeforeStr = dayBefore.strftime('%Y%m%d')
        dfSat.loc[:, 'startTime'] = dfSat.orbit.str[42:50]
        match = dfSat.loc[(dfSat.startTime == dayBeforeStr, 'orbit')].values[0]
        orbitUrl = f"{urlResorb}/{match}"
        os.system('mv ./orbits/2*/*EOF ./orbits/')

    return orbitUrl
