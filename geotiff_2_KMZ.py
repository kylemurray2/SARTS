import rasterio
import simplekml
import zipfile
import numpy as np
import matplotlib.pyplot as plt

def geotiff_to_kmz(geotiff_path, kmz_path, colormap='viridis'):
    # Step 1: Read the bounds and data from the GeoTIFF file
    with rasterio.open(geotiff_path) as src:
        bounds = src.bounds
        data = src.read(1)  # Assuming single-band image

    # Create an alpha (transparency) channel, and set it to fully opaque (255)
    alpha_channel = np.ones(data.shape, dtype=np.uint8) * 255
    
    # Set alpha to 0 (fully transparent) where data is NaN
    alpha_channel[np.isnan(data)] = 0
    
    # Rescale data to Byte (values between 0 and 255)
    data_min, data_max = np.nanmin(data), np.nanmax(data)
    norm_data = (data - data_min) / (data_max - data_min)  # Normalize data to 0-1 range
    
    # Apply colormap to normalized data
    cm = plt.get_cmap(colormap)
    colored_data = (cm(norm_data)[:, :, :3] * 255).astype(np.uint8)  # Convert colormap float values to uint8 for RGB image
    
    # Combine RGB and alpha channel to get an RGBA image
    rgba_data = np.dstack((colored_data, alpha_channel))
    
    # Step 2: Save the RGBA data to PNG
    png_path = "temp_image.png"
    with rasterio.open(png_path, 'w', driver='PNG', width=data.shape[1], height=data.shape[0], count=4, dtype=np.uint8) as dst:
        for k, channel in enumerate([rgba_data[:, :, i] for i in range(4)]):  # Split RGBA
            dst.write(channel, k + 1)


    # Step 3: Create a KML representation of the raster image
    kml = simplekml.Kml()
    ground = kml.newgroundoverlay(name="GeoTIFF as KMZ")
    ground.icon.href = png_path
    ground.draworder = 0
    ground.visibility = 1
    ground.latlonbox.north = bounds.top
    ground.latlonbox.south = bounds.bottom
    ground.latlonbox.east = bounds.right
    ground.latlonbox.west = bounds.left
    ground.latlonbox.rotation = 0
    ground.altitudemode = 'clampToSeaFloor'
    ground.gxaltitudemode ='clampToSeaFloor'
    kml.savekmz(kmz_path)




geotiff_path = "geotifs/OahuVertical_new.tif"
kmz_path = "geotifs/test.kmz"
geotiff_to_kmz(geotiff_path, kmz_path, colormap='RdBu')  # For example, 'jet' colormap