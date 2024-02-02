#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:51:21 2020


Map an IFG or other gridded data

@author: km
"""


import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter, LatitudeLocator
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
from matplotlib import pyplot as plt
import numpy as np
from cartopy.io.img_tiles import GoogleTiles
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature


class ShadedReliefESRI(GoogleTiles):
    # shaded relief
    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
               'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
               z=z, y=y, x=x)
        return url

def configure_gridlines(ax, minlon, maxlon, pad, tick_increment):
    """
    Configures the gridlines for the map.

    Parameters:
    - ax: The Axes object to apply gridlines to.
    - minlon: Minimum longitude value.
    - maxlon: Maximum longitude value.
    - pad: Padding added to longitude and latitude for extent.
    - tick_increment: Increment for longitude gridline ticks.

    Returns:
    - gl: The configured gridline object.
    """
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(np.arange(np.floor(minlon-pad), np.ceil(maxlon+pad), tick_increment))
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.top_labels = False
    gl.right_labels = False

    return gl

def plot_data(ax, lons, lats, img, contour, vmin, vmax, cm, alpha):
    """
    Plots the data on the map.

    Parameters:
    - ax: The Axes object to plot the data on.
    - lons, lats: Longitude and Latitude arrays.
    - img: The data array to plot.
    - contour: Boolean flag to choose between contourf and pcolormesh.
    - vmin, vmax: Minimum and maximum values for the color scale.
    - cm: Colormap name.
    - alpha: Alpha level for the plot.

    Returns:
    - img_handle: Handle for the plotted image, useful for creating colorbars.
    """
    cmap = plt.get_cmap(cm)

    if contour:
        img_handle = ax.contourf(lons, lats, img, levels=np.arange(vmin, vmax, 1),
                                 cmap=cmap, transform=ccrs.PlateCarree(), rasterized=True)
    else:
        img_handle = ax.pcolormesh(lons, lats, img, cmap=cmap, alpha=alpha,
                                   vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(),
                                   rasterized=True, linewidth=0, ls=":", edgecolor='face')

    return img_handle


def plot_faults(ax, fault_file):
    """
    Plots fault lines on the map.

    Parameters:
    - ax: The Axes object to plot faults on.
    - fault_file: The file path to the shapefile containing fault data.

    """

    if fault_file:
        # Read the shapefile
        reader = shpreader.Reader(fault_file)

        # Create a feature from the shapefile geometries
        shape_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=1, zorder=5, alpha=0.8)

        # Add the feature to the Axes
        ax.add_feature(shape_feature)


def add_colorbar(fig, img_handle, label, orientation='horizontal'):
    """
    Adds a colorbar to the figure.

    Parameters:
    - fig: The Figure object to which the colorbar will be added.
    - img_handle: The image handle (returned by contourf, pcolormesh, etc.) to which the colorbar is linked.
    - label: The label for the colorbar.
    - orientation: Orientation of the colorbar, 'horizontal' or 'vertical'.
    """
    cbar = fig.colorbar(img_handle, fraction=0.03, pad=0.05, orientation=orientation, label=label)
    return cbar


def mapImg(img, lons, lats, vmin, vmax, padding, zoom_level, title, background='World_Imagery', colormap='jet', figsize=(8,8), alpha=1, draw_contour=False, label='cm/yr', fault_file=None):
    """
    Plots geospatial data on a map.

    :param img: Array of image data to plot.
    :param lons: Array of longitudes.
    :param lats: Array of latitudes.
    :param vmin: Minimum value for colormap.
    :param vmax: Maximum value for colormap.
    :param padding: Padding to add around the map boundaries.
    :param zoom_level: Zoom level for the background map.
    :param title: Title of the plot.
    :param background: Background map type.
    :param colormap: Colormap for the image data.
    :param figsize: Size of the figure.
    :param alpha: Alpha value for the image plot.
    :param draw_contour: Boolean to draw contour if True.
    :param label: Label for the colorbar.
    :param fault_file: Path to a fault file to plot.
    :return: None.
    """
    # Ensure lons and lats are valid
    if len(lons) != len(lats):
        raise ValueError("Longitude and latitude arrays must be of the same length.")

    # Handling NaN values
    lons, lats = np.asarray(lons), np.asarray(lats)
    if np.isnan(lons).any() or np.isnan(lats).any():
        raise ValueError("Longitude and latitude arrays must not contain NaN values.")

    # Calculate bounds
    minlat, maxlat = np.nanmin(lats), np.nanmax(lats)
    minlon, maxlon = np.nanmin(lons), np.nanmax(lons)

    # Set up mapping tiles
    url = f'https://server.arcgisonline.com/ArcGIS/rest/services/{background}/MapServer/tile/{{z}}/{{y}}/{{x}}.jpg'
    image = cimgt.GoogleTiles(url=url)
    data_crs = image.crs

    # Initialize plot
    plt.figure(figsize=figsize)
    ax = plt.axes(projection=data_crs)
    ax.set_extent([minlon-padding, maxlon+padding, minlat-padding, maxlat+padding], crs=ccrs.PlateCarree())

    configure_gridlines(ax, minlon, maxlon, padding)
    ax.add_image(image, zoom_level)

    # Plot data
    plot_data(ax, img, lons, lats, vmin, vmax, alpha, colormap, draw_contour)
    add_colorbar(ax, label)

    # Plot faults if file is provided
    if fault_file:
        plot_faults(ax, fault_file)

    plt.title(title)
    plt.show()

    
def mapBackground(bg, minlon, maxlon, minlat, maxlat, zoomLevel, title, pad=0, scalebar=100, borders=True, fault_file=None,figsize=(8,8)):
    """
    Makes a background map to plot various data over it.
    
    Parameters:
    - bg: Background type from the ARCGIS database.
    - minlon, maxlon, minlat, maxlat: Longitude and latitude bounds.
    - zoomLevel: Zoom level for the map.
    - title: Title of the map.
    - pad: Padding around the bounds.
    - scalebar: Length of the scale bar.
    - borders: Whether to plot country borders.
    - fault_file: path/file.shp shp file for faults (optional)
    - figsize: fig size (x,y)
    Background options:
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
    """
    url = 'https://server.arcgisonline.com/ArcGIS/rest/services/' + bg + '/MapServer/tile/{z}/{y}/{x}.jpg'
    image = cimgt.GoogleTiles(url=url)
    data_crs = image.crs

    plt.figure(figsize=figsize)
    ax = plt.axes(projection=data_crs)
    ax.set_extent([minlon-pad, maxlon+pad, minlat-pad, maxlat+pad], crs=ccrs.PlateCarree())

    if borders:
        ax.add_feature(cfeature.BORDERS, linewidth=.5, color='grey')
        ax.coastlines(linewidth=.5)


    # Calculate tick increment for gridlines
    lon_range = (pad+maxlon) - (minlon-pad)
    lat_range = (pad+maxlat) - (minlat-pad)
    range_min = min(lon_range, lat_range)
    tick_increment = round(range_min / 4, 1)

    # Configure gridlines
    configure_gridlines(ax, minlon, maxlon, pad, tick_increment)

    # Add background image
    ax.add_image(image, zoomLevel, zorder=1)

    # Plot faults if required
    if fault_file:
        plot_faults(ax, fault_file)  # You will need to provide the appropriate fault file

    if scalebar:  
        scale_bar(ax, (0.1,.1), scalebar)


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
        return geodesic.inverse(a_phys, b_phys)[0][0]

    end = _upper_bound(start, direction, distance, dist_func)

    return _distance_along_line(start, end, distance, dist_func, tol)


def scale_bar(ax, location, length, metres_per_unit=1000, unit_name='km',
              tol=0.01, angle=0, color='black', linewidth=3, text_offset=0.005,
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
    end = _point_along_line(ax, location, length_metres, angle=angle_rad,
                            tol=tol)

    # Coordinates are currently in axes coordinates, so use transAxes to
    # put into data coordinates. *zip(a, b) produces a list of x-coords,
    # then a list of y-coords.
    ax.plot(*zip(location, end), transform=ax.transAxes, **plot_kwargs)

    # Push text away from bar in the perpendicular direction.
    midpoint = (location + end) / 2
    offset = text_offset * np.array([-np.sin(angle_rad), np.cos(angle_rad)])
    text_location = midpoint + offset

    # 'rotation' keyword argument is in text_kwargs.
    ax.text(*text_location, f"{length} {unit_name}", rotation_mode='anchor',
            transform=ax.transAxes,fontsize=8, **text_kwargs)