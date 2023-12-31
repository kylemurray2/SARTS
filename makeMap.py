#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:51:21 2020


Map an IFG or other gridded data

@author: km
"""


from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.img_tiles as cimgt
from cartopy import config
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
from matplotlib import pyplot as plt
import numpy as np
from cartopy.io.img_tiles import GoogleTiles


class ShadedReliefESRI(GoogleTiles):
    # shaded relief
    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
               'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
               z=z, y=y, x=x)
        return url


def mapImg(img, lons, lats, vmin, vmax, pad,zoom, title, bg='World_Imagery', cm='jet', plotFaults= False,alpha=1,contour=False):
    
    minlat=lats.min()
    maxlat=lats.max()
    minlon=lons.min()
    maxlon=lons.max()
    url = 'https://server.arcgisonline.com/ArcGIS/rest/services/' + bg + '/MapServer/tile/{z}/{y}/{x}.jpg'
    image = cimgt.GoogleTiles(url=url)
    data_crs = image.crs #ShadedReliefESRI().crs#ccrs.PlateCarree()
    fig =  plt.figure(figsize=(8,8))
    ax = plt.axes(projection=data_crs)
    ax.set_extent([minlon-pad, maxlon+pad, minlat-pad, maxlat+pad], crs=ccrs.PlateCarree())
    cmap = plt.get_cmap(cm)

    
    lon_range = (pad+maxlon) - (minlon-pad)
    lat_range = (pad+maxlat) - (minlat-pad)
    rangeMin = np.min(np.array([lon_range,lat_range]))
    tick_increment = round(rangeMin/8,2)
    
    import matplotlib.ticker as mticker
    from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
    
    # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    #               linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    # gl.xlocator = mticker.FixedLocator(np.arange(np.floor(minlon-pad),np.ceil(maxlon+pad),tick_increment))
    # gl.ylocator = LatitudeLocator()
    # gl.xformatter = LongitudeFormatter()
    # gl.yformatter = LatitudeFormatter()
    # gl.ylabel_style = {'size': 8, 'color': 'black'}
    # gl.xlabel_style = {'size': 8, 'color': 'black'}
    # gl.top_labels = False
    # gl.right_labels = False
    
    ax.add_image(image, zoom) #zoom level
    
    if contour:
        img_handle = plt.contourf(lons, lats,img,levels=np.arange(vmin,vmax,1),cmap=cmap,transform=ccrs.PlateCarree(),rasterized = True)
    else:
        img_handle = plt.pcolormesh(lons, lats, img,cmap=cmap, alpha=alpha,vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree(),rasterized = True,linewidth=0,ls=":", edgecolor='face')


    plt.colorbar(img_handle,fraction=0.03, pad=0.05,orientation='horizontal',label='mm/yr')
    
    if plotFaults:
        # Plot faults
        import cartopy.io.shapereader as shpreader
        from cartopy.feature import ShapelyFeature
        # reader = shpreader.Reader("/d/MapData/gem-global-active-faults/shapefile/gem_active_faults.shp")
        reader = shpreader.Reader("/d/MapData/EARS/kivu.shp")

        shape_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), edgecolor='r', facecolor='none',linewidth=1,zorder=5,alpha=0.8)
        ax.add_feature(shape_feature)
    
    plt.title(title)
    plt.show()
    
    
def mapBackground(bg, minlon, maxlon, minlat, maxlat, zoomLevel, title, pad=0, scalebar=100, borders=True,plotFaults=True):
  
    '''
    Makes a background map that you can then plot stuff over (footprints, scatterplot, etc.)
    bg: background type from the ARCGIS database choose from the following:
        NatGeo_world_Map
        USA_Topo_Maps
        World_Imagery
        World_Physical_Map
        World_Shaded_Relief 
        World_Street_Map 
        World_Terrain_Base 
        World_Topo_Map 
        Specialty/DeLorme_World_Base_Map
        Specialty/World_Navigation_Charts
        Canvas/World_Dark_Gray_Base 
        Canvas/World_Dark_Gray_Reference
        Canvas/World_Light_Gray_Base 
        Canvas/World_Light_Gray_Reference
        Elevation/World_Hillshade_Dark 
        Elevation/World_Hillshade 
        Ocean/World_Ocean_Base 
        Ocean/World_Ocean_Reference 
        Polar/Antarctic_Imagery
        Polar/Arctic_Imagery
        Polar/Arctic_Ocean_Base
        Polar/Arctic_Ocean_Reference
        Reference/World_Boundaries_and_Places_Alternate
        Reference/World_Boundaries_and_Places
        Reference/World_Reference_Overlay
        Reference/World_Transportation 
        WorldElevation3D/Terrain3D
        WorldElevation3D/TopoBathy3D
    '''
    
    url = 'https://server.arcgisonline.com/ArcGIS/rest/services/' + bg + '/MapServer/tile/{z}/{y}/{x}.jpg'
    image = cimgt.GoogleTiles(url=url)
    data_crs = image.crs #ShadedReliefESRI().crs#ccrs.PlateCarree()
    fig =  plt.figure(figsize=(6,6))
    ax = plt.axes(projection=data_crs)
    ax.set_extent([minlon-pad, maxlon+pad, minlat-pad, maxlat+pad], crs=ccrs.PlateCarree())

    if borders:
        ax.add_feature(cfeature.BORDERS,linewidth=1,color='red')
    # ax.add_feature(cfeature.OCEAN)
    # ax.add_feature(cfeature.LAKES)
    # ax.add_feature(cfeature.RIVERS)
    
    
    lon_range = (pad+maxlon) - (minlon-pad)
    lat_range = (pad+maxlat) - (minlat-pad)
    rangeMin = np.min(np.array([lon_range,lat_range]))
    tick_increment = round(rangeMin/4,1)
    
    import matplotlib.ticker as mticker
    from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(np.arange(np.floor(minlon-pad),np.ceil(maxlon+pad),tick_increment))
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.top_labels = False
    gl.right_labels = False
    
    ax.add_image(image, zoomLevel,zorder=1) #zoom level
    
    if plotFaults:
        # Plot faults
        import cartopy.io.shapereader as shpreader
        from cartopy.feature import ShapelyFeature
        if minlon >0:
            # reader = shpreader.Reader("/d/MapData/gem-global-active-faults/shapefile/gem_active_faults.shp")
            reader = shpreader.Reader("/d/MapData/EARS/kivu.shp")

        else:
            # reader2 = shpreader.Reader("/d/faults/CFM/traces/shp/CFM5.3_traces.shp")
            reader = shpreader.Reader("/d/faults/Shapefile/QFaults.shp")
        
        shape_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), edgecolor='r', facecolor='none',linewidth=.5,zorder=5,alpha=0.8)
        ax.add_feature(shape_feature)
        # shape_feature2 = ShapelyFeature(reader2.geometries(), ccrs.PlateCarree(), edgecolor='b', facecolor='none',linewidth=.5,zorder=5,alpha=0.8)
        # ax.add_feature(shape_feature2)
    
    # scale_bar(ax, location, length, metres_per_unit=1000, unit_name='km',
    #                tol=0.01, angle=0, color='black', linewidth=5, text_offset=0.01,
    #                ha='center', va='bottom', plot_kwargs=None, text_kwargs=None)
    if scalebar:
        scale_bar(ax, (.1,.1), scalebar,linewidth=0.5)
    

    plt.title(title)
    plt.show()
    return data_crs






def _axes_to_lonlat(ax, coords):
    """(lon, lat) from axes coordinates."""
    display = ax.transAxes.transform(coords)
    data = ax.transData.inverted().transform(display)
    lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)

    return lonlat


def _upper_bound(start, direction, distance, dist_func):
    """A point farther than distance from start, in the given direction.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        direction  Nonzero (2, 1)-shaped array, a direction vector.
        distance:  Positive distance to go past.
        dist_func: A two-argument function which returns distance.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    if distance <= 0:
        raise ValueError(f"Minimum distance is not positive: {distance}")

    if np.linalg.norm(direction) == 0:
        raise ValueError("Direction vector must not be zero.")

    # Exponential search until the distance between start and end is
    # greater than the given limit.
    length = 0.1
    end = start + length * direction

    while dist_func(start, end) < distance:
        length *= 2
        end = start + length * direction

    return end


def _distance_along_line(start, end, distance, dist_func, tol):
    """Point at a distance from start on the segment  from start to end.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        end:       Outer bound on point's location.
        distance:  Positive distance to travel.
        dist_func: Two-argument function which returns distance.
        tol:       Relative error in distance to allow.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    initial_distance = dist_func(start, end)
    if initial_distance < distance:
        raise ValueError(f"End is closer to start ({initial_distance}) than "
                         f"given distance ({distance}).")

    if tol <= 0:
        raise ValueError(f"Tolerance is not positive: {tol}")

    # Binary search for a point at the given distance.
    left = start
    right = end

    while not np.isclose(dist_func(start, right), distance, rtol=tol):
        midpoint = (left + right) / 2

        # If midpoint is too close, search in second half.
        if dist_func(start, midpoint) < distance:
            left = midpoint
        # Otherwise the midpoint is too far, so search in first half.
        else:
            right = midpoint

    return right


def _point_along_line(ax, start, distance, angle=0, tol=0.01):
    """Point at a given distance from start at a given angle.

    Args:
        ax:       CartoPy axes.
        start:    Starting point for the line in axes coordinates.
        distance: Positive physical distance to travel.
        angle:    Anti-clockwise angle for the bar, in radians. Default: 0
        tol:      Relative error in distance to allow. Default: 0.01

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    # Direction vector of the line in axes coordinates.
    direction = np.array([np.cos(angle), np.sin(angle)])

    geodesic = cgeo.Geodesic()

    # Physical distance between points.
    def dist_func(a_axes, b_axes):
        a_phys = _axes_to_lonlat(ax, a_axes)
        b_phys = _axes_to_lonlat(ax, b_axes)

        # Geodesic().inverse returns a NumPy MemoryView like [[distance,
        # start azimuth, end azimuth]].
        return geodesic.inverse(a_phys, b_phys).base[0, 0]

    end = _upper_bound(start, direction, distance, dist_func)

    return _distance_along_line(start, end, distance, dist_func, tol)


def scale_bar(ax, location, length, metres_per_unit=1000, unit_name='km',
              tol=0.01, angle=0, color='black', linewidth=5, text_offset=0.01,
              ha='center', va='bottom', plot_kwargs=None, text_kwargs=None,
              **kwargs):
    """Add a scale bar to CartoPy axes.

    For angles between 0 and 90 the text and line may be plotted at
    slightly different angles for unknown reasons. To work around this,
    override the 'rotation' keyword argument with text_kwargs.

    Args:
        ax:              CartoPy axes.
        location:        Position of left-side of bar in axes coordinates.
        length:          Geodesic length of the scale bar.
        metres_per_unit: Number of metres in the given unit. Default: 1000
        unit_name:       Name of the given unit. Default: 'km'
        tol:             Allowed relative error in length of bar. Default: 0.01
        angle:           Anti-clockwise rotation of the bar.
        color:           Color of the bar and text. Default: 'black'
        linewidth:       Same argument as for plot.
        text_offset:     Perpendicular offset for text in axes coordinates.
                         Default: 0.005
        ha:              Horizontal alignment. Default: 'center'
        va:              Vertical alignment. Default: 'bottom'
        **plot_kwargs:   Keyword arguments for plot, overridden by **kwargs.
        **text_kwargs:   Keyword arguments for text, overridden by **kwargs.
        **kwargs:        Keyword arguments for both plot and text.
    """
    # Setup kwargs, update plot_kwargs and text_kwargs.
    if plot_kwargs is None:
        plot_kwargs = {}
    if text_kwargs is None:
        text_kwargs = {}

    plot_kwargs = {'linewidth': linewidth, 'color': color, **plot_kwargs,
                   **kwargs}
    text_kwargs = {'ha': ha, 'va': va, 'rotation': angle, 'color': color,
                   **text_kwargs, **kwargs}

    # Convert all units and types.
    location = np.asarray(location)  # For vector addition.
    length_metres = length * metres_per_unit
    angle_rad = angle * np.pi / 180

    # End-point of bar.
    end = _point_along_line(ax, location, length_metres, angle=angle_rad,tol=tol)

    from matplotlib import patheffects
    buffer = [patheffects.withStroke(linewidth=1, foreground="w")]
    # Coordinates are currently in axes coordinates, so use transAxes to
    # put into data coordinates. *zip(a, b) produces a list of x-coords,
    # then a list of y-coords.
    ax.plot(*zip(location, end), transform=ax.transAxes, linewidth=linewidth, color='black',  alpha=.6)

    # Push text away from bar in the perpendicular direction.
    midpoint = (location + end) / 2
    offset = text_offset * np.array([-np.sin(angle_rad), np.cos(angle_rad)])
    text_location = midpoint + offset

    # 'rotation' keyword argument is in text_kwargs.
    ax.text(*text_location, f"{length} {unit_name}", rotation_mode='anchor',
            transform=ax.transAxes, **text_kwargs,  path_effects=buffer)

