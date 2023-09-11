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




def fetch_orbit(granuleName, base_url):
    sat = granuleName[:3]
    date = granuleName[17:25]
    print(f"Retrieving precise orbit URL for {sat}, {date}")

    try:
        r = requests.get(base_url, stream=True)
        r.raise_for_status()
    except requests.RequestException as e:
        print(f"Failed to fetch data: {e}")
        return None

    orbits = html.fromstring(r.content).xpath('//a/@href')
    return find_matching_orbit(orbits, sat, date)

def find_matching_orbit(orbits, sat, date):
    df = pd.DataFrame({'orbit': orbits})
    filtered_orbits = df[df['orbit'].str.startswith(sat)].copy()
    previous_day = pd.to_datetime(date) - timedelta(days=1)
    previous_day_str = previous_day.strftime('%Y%m%d')

    filtered_orbits['startTime'] = filtered_orbits['orbit'].str[42:50]
    matched_orbit = filtered_orbits[filtered_orbits['startTime'] == previous_day_str]['orbit'].values

    if matched_orbit:
        return matched_orbit[0]
    else:
        print(f"No matching orbit found for {sat} on {previous_day_str}")
        return None

def get_orbit_url(granuleName):
    """Retrieve precise orbit file for a specific Sentinel-1 granule."""

    precise_url = 'https://s1qc.asf.alaska.edu/aux_poeorb'
    resorb_url = 'https://s1qc.asf.alaska.edu/aux_resorb'

    matched_orbit = fetch_orbit(granuleName, precise_url)

    if matched_orbit:
        return f"{precise_url}/{matched_orbit}"
    else:
        print('Using resorb for this date (it is probably too recent).')
        matched_orbit = fetch_orbit(granuleName, resorb_url)

        if matched_orbit:
            os.system('mv ./orbits/2*/*EOF ./orbits/')
            return f"{resorb_url}/{matched_orbit}"

    return None






# def get_orbit_url(granuleName):
#     """Retrieve precise orbit file for a specific Sentinel-1 granule.
#     Precise orbits available ~3 weeks after aquisition.
#     Parameters
#     ----------
#     granuleName : str
#         ASF granule name, e.g.:
#         S1B_IW_SLC__1SDV_20171117T015310_20171117T015337_008315_00EB6C_40CA
#     url : str
#         website with simple list of orbit file links
#     Returns
#     -------
#     orbitUrl :  str
#         url pointing to matched orbit file
#     """

#     urlPrecise='https://s1qc.asf.alaska.edu/aux_poeorb'
#     urlResorb='https://s1qc.asf.alaska.edu/aux_resorb'

#     sat = granuleName[:3]
#     date = granuleName[17:25]
#     print(f"retrieving precise orbit URL for {sat}, {date}")
#     try:
#         r = requests.get(urlPrecise,stream=True)
#         webpage = html.fromstring(r.content)
#         orbits = webpage.xpath('//a/@href')
#         df = pd.DataFrame(dict(orbit=orbits))
#         dfSat = df[df.orbit.str.startswith(sat)].copy()
#         dayBefore = pd.to_datetime(date) - pd.to_timedelta(1, unit='d')
#         dayBeforeStr = dayBefore.strftime('%Y%m%d')
#         dfSat.loc[:, 'startTime'] = dfSat.orbit.str[42:50]
#         match = dfSat.loc[(dfSat.startTime == dayBeforeStr, 'orbit')].values[0]
#         orbitUrl = f"{urlPrecise}/{match}"
#     except:
#         print('using resorb for this date (it is probably too recent)')
#         r = requests.get(urlResorb)
#         webpage = html.fromstring(r.content)
#         orbits = webpage.xpath('//a/@href')
#         df = pd.DataFrame(dict(orbit=orbits))
#         dfSat = df[df.orbit.str.startswith(sat)].copy()
#         dayBefore = pd.to_datetime(date) - pd.to_timedelta(1, unit='d')
#         dayBeforeStr = dayBefore.strftime('%Y%m%d')
#         dfSat.loc[:, 'startTime'] = dfSat.orbit.str[42:50]
#         match = dfSat.loc[(dfSat.startTime == dayBeforeStr, 'orbit')].values[0]
#         orbitUrl = f"{urlResorb}/{match}"
#         os.system('mv ./orbits/2*/*EOF ./orbits/')

#     return orbitUrl
