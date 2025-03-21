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

def getGran(path, start, end, sat, bounds, poly, processing_level='BURST'):
    """
    Query ASF for granules.
    
    Parameters
    ----------
    path : str
        Relative orbit path
    start : str
        Start date in format 'YYYY-MM-DDThh:mm:ssZ'
    end : str
        End date in format 'YYYY-MM-DDThh:mm:ssZ'
    sat : str
        Satellite platform (e.g., 'S1' for Sentinel-1)
    bounds : str
        Bounds in format 'S,N,W,E'
    poly : str
        WKT polygon string
    processing_level : str
        Processing level (default: 'BURST')
        
    Returns
    -------
    tuple
        (URLs, granule names, acquisition dates, response object)
    """
    fmt = 'CSV'
    flightDirection = None
    baseurl = 'https://api.daac.asf.alaska.edu/services/search/param'
    
    # Set up basic query parameters
    data = dict(
        platform=sat,  # Use 'S1' or 'Sentinel-1' for Sentinel-1 data
        output=fmt
    )

    # Add the processing level
    data['processingLevel'] = processing_level

    # Handle special case for ALOS
    if sat == 'ALOS':
        data['beamMode'] = 'FBS,FBD'
        data['processingLevel'] = 'L1.0'
    
    # Set up spatial search parameters
    if poly:
        print('Using polygon')
        data['intersectsWith'] = poly
    elif bounds:
        print('Using bounds. Usually it is better to use a polygon instead.')
        miny, maxy, minx, maxx = bounds.split(sep=',')
        roi = shapely.geometry.box(float(minx), float(miny), float(maxx), float(maxy))
        polygonWKT = roi.wkt
        data['intersectsWith'] = polygonWKT
    else:
        print('Need to specify the polygon search area')
        sys.exit(1)

    # Add optional parameters if provided
    if path:
        data['relativeOrbit'] = path
    if start:
        data['start'] = start
    if end:
        data['end'] = end
    if flightDirection:
        data['flightDirection'] = flightDirection

    print('Searching for data with the following parameters:')
    for key, value in data.items():
        if key != 'intersectsWith':  # Don't print the full polygon WKT to keep output clean
            print(f"  {key}: {value}")
        else:
            print(f"  {key}: [WKT polygon]")

    # Make the API request
    print('Sending request to ASF API...')
    r = requests.get(baseurl, params=data, timeout=100)
    print(f"Request URL: {r.url}")
    
    # Save raw response to file
    with open('out.csv', 'w') as j:
        j.write(r.text)
    
    # Check if the response is empty or contains an error
    if len(r.text.strip()) == 0 or "Not Found" in r.text:
        print("Error: Empty or error response received from ASF API")
        return ([], [], [], r)
    
    try:
        # Parse the CSV response
        df = pd.read_csv('out.csv')
        
        if df.empty:
            print("No data found matching your search criteria.")
            return ([], [], [], r)
        
        # Extract relevant columns
        slcUrls = df['URL']
        
        if 'Size (MB)' in df.columns:
            sizesMB = df['Size (MB)']
            print(f'Found {round((sum(sizesMB)/1000),2)} GB of data')
        else:
            print('Size information not available in the response')
        
        if 'Acquisition Date' in df.columns:
            dates = df['Acquisition Date']
        else:
            dates = [None] * len(slcUrls)
            print('Acquisition date information not available in the response')
        
        if 'Granule Name' in df.columns:
            gran = df['Granule Name']
        else:
            gran = [None] * len(slcUrls)
            print('Granule name information not available in the response')
        
        print(f"Found {len(slcUrls)} results")
        
        # Display some example results if available
        if len(slcUrls) > 0:
            print("\nExample results:")
            for i in range(min(3, len(slcUrls))):
                print(f"  {i+1}. {gran.iloc[i] if gran is not None else 'N/A'} - {dates.iloc[i] if dates is not None else 'N/A'}")
            print("...")
        
        return (slcUrls, gran, dates, r)
    
    except Exception as e:
        print(f"Error parsing response: {e}")
        print("Response content:")
        print(r.text[:500] + "..." if len(r.text) > 500 else r.text)
        return ([], [], [], r)


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


def get_orbit_url(granuleNames):
    """Retrieve precise orbit file for a specific Sentinel-1 granule or date.
    Precise orbits available ~3 weeks after acquisition.
    
    Parameters
    ----------
    granuleNames : list
        List of granule names or satellite+date strings (e.g., S1A20200229)
        Standard SLC granule: S1B_IW_SLC__1SDV_20171117T015310_20171117T015337_008315_00EB6C_40CA
        Burst format: S1_322422_IW3_20200229T010154_VH_C8E8-BURST
        Or direct format: S1A20200229
        
    Returns
    -------
    orbitUrl : list
        List of URLs pointing to matched orbit files
    """
    urlPrecise = 'https://s1qc.asf.alaska.edu/aux_poeorb'
    urlResorb = 'https://s1qc.asf.alaska.edu/aux_resorb'
    
    print("Fetching available orbit files...")
    # Get precise orbit file list
    r = get_with_retry(urlPrecise)
    webpage = html.fromstring(r.content)
    orbits = webpage.xpath('//a/@href')
    df = pd.DataFrame(dict(orbit=orbits))

    # Get restituted orbit file list
    rr = get_with_retry(urlResorb)
    webpager = html.fromstring(rr.content)
    orbitsr = webpager.xpath('//a/@href')
    dfr = pd.DataFrame(dict(orbit=orbitsr))

    orbitUrls = []
    
    for granuleName in granuleNames:
        print(f"Processing: {granuleName}")
        
        # Case 1: Direct satellite+date format (e.g., S1A20200229)
        if len(granuleName) >= 11 and granuleName.startswith('S1'):
            # Already in the format we need
            sat = granuleName[:3]  # S1A or S1B
            date_str = granuleName[3:11]  # YYYYMMDD
        
        # Case 2: Check if this is a burst format granule with date
        elif "-BURST" in granuleName:
            parts = granuleName.split('_')
            # Find the part that looks like a date (starts with year)
            date_part = next((part for part in parts if len(part) >= 8 and part[:4].isdigit() and 'T' in part), None)
            
            if date_part:
                # Format is like 20200229T010154
                date_str = date_part[:8]  # Extract YYYYMMDD part
                
                # For bursts, S1A/S1B is not always clear, so check for indicators or default to S1A
                if "1A_" in granuleName or "1A-" in granuleName:
                    sat = "S1A"
                elif "1B_" in granuleName or "1B-" in granuleName:
                    sat = "S1B"
                else:
                    # Default to S1A if not specified (can be refined if needed)
                    sat = "S1A"
            else:
                print(f"  Cannot parse date from burst format: {granuleName}")
                orbitUrls.append(None)
                continue
        
        # Case 3: Standard SLC granule format
        elif granuleName.startswith('S1') and '_' in granuleName:
            try:
                sat = granuleName[:3]  # Extract S1A or S1B
                date_parts = [p for p in granuleName.split('_') if len(p) >= 8 and p[0:4].isdigit() and 'T' in p]
                if date_parts:
                    date_str = date_parts[0][:8]  # Extract YYYYMMDD part from first date component
                else:
                    print(f"  Cannot parse date from standard format: {granuleName}")
                    orbitUrls.append(None)
                    continue
            except:
                print(f"  Error parsing standard granule name: {granuleName}")
                orbitUrls.append(None)
                continue
        
        else:
            print(f"  Unrecognized granule format, cannot parse: {granuleName}")
            orbitUrls.append(None)
            continue
        
        print(f"  Using satellite: {sat}, date: {date_str}")
        
        try:
            # Try to find precise orbit file first
            dfSat = df[df.orbit.str.startswith(sat)].copy()
            if dfSat.empty:
                print(f"  No precise orbit files found for {sat}")
                raise ValueError("No precise orbit files for this satellite")
                
            dayBefore = pd.to_datetime(date_str) - pd.to_timedelta(1, unit='d')
            dayBeforeStr = dayBefore.strftime('%Y%m%d')
            
            dfSat.loc[:, 'startTime'] = dfSat.orbit.str[42:50]
            matches = dfSat.loc[dfSat.startTime == dayBeforeStr, 'orbit'].values
            
            if len(matches) > 0:
                match = matches[0]
                orbitUrl = f"{urlPrecise}/{match}"
                print(f"  Found precise orbit file: {match}")
            else:
                raise ValueError(f"No precise orbit file matches date {dayBeforeStr}")
                
        except Exception as e:
            try:
                # Try restituted orbit file as fallback
                print(f"  Trying restituted orbit for {sat} {date_str} (reason: {str(e)})")
                dfSat = dfr[dfr.orbit.str.startswith(sat)].copy()
                
                if dfSat.empty:
                    print(f"  No restituted orbit files found for {sat}")
                    orbitUrls.append(None)
                    continue
                    
                dayBefore = pd.to_datetime(date_str) - pd.to_timedelta(1, unit='d')
                dayBeforeStr = dayBefore.strftime('%Y%m%d')
                
                dfSat.loc[:, 'startTime'] = dfSat.orbit.str[42:50]
                matches = dfSat.loc[dfSat.startTime == dayBeforeStr, 'orbit'].values
                
                if len(matches) > 0:
                    match = matches[0]
                    orbitUrl = f"{urlResorb}/{match}"
                    print(f"  Found restituted orbit file: {match}")
                else:
                    print(f"  No orbit file found for {sat} {date_str}")
                    orbitUrls.append(None)
                    continue
            except Exception as nested_e:
                print(f"  Failed to find any orbit file for {sat} {date_str}: {nested_e}")
                orbitUrls.append(None)
                continue
        
        orbitUrls.append(orbitUrl)
    
    # Report statistics
    valid_orbits = sum(1 for url in orbitUrls if url is not None)
    print(f"Found {valid_orbits} orbit files out of {len(granuleNames)} requested granules")
    
    return orbitUrls


def save_query_results(urls, granules, dates, output_file='query_results.csv'):
    """
    Save the query results to a CSV file.
    
    Parameters
    ----------
    urls : list
        List of URLs for burst files
    granules : list
        List of granule names
    dates : list
        List of acquisition dates
    output_file : str
        Output CSV filename
    """
    df = pd.DataFrame({
        'URL': urls,
        'Granule Name': granules,
        'Acquisition Date': dates
    })
    df.to_csv(output_file, index=False)
    print(f"Query results saved to {output_file}")
    return df


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Query for Sentinel-1 burst data from ASF')
    parser.add_argument('--path', type=str, help='Relative orbit path')
    parser.add_argument('--start', type=str, help='Start date (format: YYYY-MM-DDThh:mm:ssZ)')
    parser.add_argument('--end', type=str, help='End date (format: YYYY-MM-DDThh:mm:ssZ)')
    parser.add_argument('--sat', type=str, default='Sentinel-1', help='Satellite (default: Sentinel-1)')
    parser.add_argument('--bounds', type=str, help='Bounds in format: S,N,W,E')
    parser.add_argument('--polygon', type=str, help='WKT polygon')
    parser.add_argument('--output', type=str, default='query_results.csv', 
                        help='Output file for query results (default: query_results.csv)')
    parser.add_argument('--orbit_info', action='store_true', 
                        help='Get orbit files information for the granules')
    parser.add_argument('--processing_level', type=str, default='BURST',
                        help='Processing level (default: BURST)')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug output')
    
    args = parser.parse_args()
    
    # Check required arguments
    if not args.polygon and not args.bounds:
        print("Error: Either --polygon or --bounds must be specified")
        sys.exit(1)
        
    if not args.start or not args.end:
        print("Error: Both --start and --end dates must be specified")
        sys.exit(1)
    
    # Query for burst data
    urls, granules, dates, response = getGran(
        path=args.path,
        start=args.start,
        end=args.end,
        sat=args.sat,
        bounds=args.bounds,
        poly=args.polygon,
        processing_level=args.processing_level
    )
    
    # Save the query results if any were found
    if len(urls) > 0:
        df = save_query_results(urls, granules, dates, args.output)
        
        # Get orbit information if requested
        if args.orbit_info and len(granules) > 0:
            print("\nGetting orbit information for granules...")
            
            # Extract unique dates from the granules for debugging
            if args.debug:
                unique_dates = set()
                for granule in granules:
                    try:
                        parts = granule.split('_')
                        date_part = next((part for part in parts if len(part) >= 8 and part[:4].isdigit()), None)
                        if date_part:
                            unique_dates.add(date_part[:8])
                    except:
                        pass
                print(f"Found unique dates in granules: {unique_dates}")
            
            orbit_urls = get_orbit_url(granules)
            
            # Add orbit URLs to the dataframe and save
            df['Orbit URL'] = orbit_urls
            df.to_csv(args.output, index=False)
            print(f"Added orbit information to {args.output}")
        
        print(f"\nQuery complete. Results saved to {args.output}")
    else:
        print("\nNo results found. Check your search parameters and try again.")
        print("\nTips for effective searching:")
        print("1. For Sentinel-1 bursts, use --sat 'Sentinel-1'")
        print("2. Make sure your date range is not too narrow")
        print("3. Verify that your polygon or bounds intersect with Sentinel-1 coverage")
        print("4. Check that the orbit path is correct for your area of interest")
