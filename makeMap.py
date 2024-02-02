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
from cartopy import config
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
    fig = plt.figure(figsize=figsize)
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

    
def mapBackground(bg, minlon, maxlon, minlat, maxlat, zoomLevel, title, pad=0, scalebar=100, borders=True, fault_file=None):
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

    fig = plt.figure(figsize=(6,6))
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
    print('hi')
    print(ax.transAxes)
    # Plot faults if required
    if fault_file:
        plot_faults(ax, fault_file)  # You will need to provide the appropriate fault file

    # Add scale bar (optional)
    # scale_bar(ax, location, length, metres_per_unit=1000, unit_name='km',
    #                tol=0.01, angle=0, color='black', linewidth=5, text_offset=0.01,
    #                ha='center', va='bottom', plot_kwargs=None, text_kwargs=None)
    # if scalebar:
        # Implement a scale_bar function or use an existing library function to add a scale bar
        # scale_bar(ax, (.1, .1), scalebar, linewidth=0.5)
    
    add_scale_bar(ax, data_crs,1000)  # Adding a 100 km scale bar

    plt.title(title)
    plt.show()
    return data_crs



def _axes_to_lonlat(ax, coords):
    """(lon, lat) from axes coordinates."""
    display = ax.transAxes.transform(coords)
    data = ax.transData.inverted().transform(display)
    lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)

    return lonlat



def add_scale_bar(ax,crs, scale_length_km, location=(0.05, 0.05)):
    """
    Adds a scale bar to a map.

    Args:
    - ax: Matplotlib Axes object to add the scale bar to.
    - scale_length_km: Length of the scale bar in kilometers.
    
    """
    import matplotlib.patches as mpatches

    # Convert the scale length from kilometers to map units (assuming the map is in meters)
    scale_length_map_units = scale_length_km * 1000  # 1 km = 1000 m

    # Get the axes bounds and calculate scale bar size and position
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    xtotallength = xlim[1]-xlim[0]
    ytotallength = ylim[1]-ylim[0]

    
    
    scale_bar_length = (scale_length_map_units / (xlim[1] - xlim[0]))
    
    scale_bar_x = xlim[0] + .05*(xtotallength) + (.5*scale_bar_length*xtotallength)
    scale_bar_y =  ylim[0] + .02*(ytotallength) + (.1*scale_bar_length*ytotallength)


    # Create a rectangle patch for scale bar
    scale_bar = mpatches.Rectangle((scale_bar_x, scale_bar_y), scale_bar_length, 0.01, 
                                    transform=crs, color='white', ec='black', lw=2)

    # Add the scale bar to the axes
    ax.add_patch(scale_bar)

    # Add text label for the scale bar
    plt.text(scale_bar_x + scale_bar_length / 2, scale_bar_y + 0.01, 
              f'{scale_length_km} km', transform=crs, ha='center', 
              va='bottom', backgroundcolor='white',zorder=14)

