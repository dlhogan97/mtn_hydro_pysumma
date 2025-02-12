# %%
import geopandas as gpd 
from pathlib import Path
import xarray as xr
import rioxarray as rxr

def clip_geopotential(basin):
    # %%
    root_path = Path(f"/storage/dlhogan/summa_modeling_data/")

    # %%
    # open the geopotential data and take the temporal average
    gph_ds = xr.open_dataset(root_path / "ERA5_surface_geopotential.nc").mean(dim="valid_time")['z']
    # add the crs to the geopotential data
    gph_ds = gph_ds.rio.write_crs("epsg:4326")

    # %%
    # open the basin shapefile
    basin_gdf = gpd.read_file(root_path / f"domain_{basin}" / "shapefiles" / "catchment" / f"{basin}.shp")

    # %%
    # clip box the geopotential data to the basin
    gph_ds_clipped = gph_ds.rio.clip_box(basin_gdf.total_bounds[0]-0.25, basin_gdf.total_bounds[1]-0.25, basin_gdf.total_bounds[2]+0.25, basin_gdf.total_bounds[3]+0.25)

    # %%
    gravity = 9.80665 # m/s^2

    # calculate the geopotential height
    gph_height = gph_ds_clipped / gravity

    # add units
    gph_height.attrs['units'] = 'm'

    # %%
    # save the geopotential height data
    gph_height.to_netcdf(root_path / f"domain_{basin}" / "forcing" / "0_geopotential" / f"ERA5_surface_geopotential_height_{basin}.nc")
    return

# %%
basin = "TuolumneRiver"
clip_geopotential(basin)

basin = "EastRiver"
clip_geopotential(basin)