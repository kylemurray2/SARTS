import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.geometry import Polygon

fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.set_extent([-180, 180, -90, 90])  # Global extent
ax.coastlines()
ax.gridlines(draw_labels=True)

coords = []

def onclick(event):
    # Check if zoom tool or any other non-default tool is active
    if fig.canvas.toolbar.mode != '':
        return

    ix, iy = event.xdata, event.ydata

    # Check distance to the start point
    if coords and abs(ix - coords[0][0]) < 1 and abs(iy - coords[0][1]) < 1:
        # Close the polygon if click is on/near the starting point
        coords.append(coords[0])
        ax.plot([coords[-2][0], coords[-1][0]], [coords[-2][1], coords[-1][1]], 'ro-')
        fig.canvas.mpl_disconnect(cid)
        polygon = Polygon(coords)
        global wkt_output
        wkt_output = polygon.wkt
        print(wkt_output)
    else:
        coords.append((ix, iy))
        if len(coords) == 1:
            # First click just plots the point
            ax.plot(ix, iy, 'ro')
        else:
            # Subsequent clicks plot lines
            ax.plot([coords[-2][0], coords[-1][0]], [coords[-2][1], coords[-1][1]], 'ro-')
    

    # Explicitly redraw the canvas after the update
    fig.canvas.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()
