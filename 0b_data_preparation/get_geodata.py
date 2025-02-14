# %%
# Script to retrieve geospatial data for a given basin, including:
# - Basin geometry
# - Topographic data (DEM, slope, aspect)
# - Land cover data
# - Streamflow observations
# - Flowlines (main stem and tributaries)
# - Subcatchments
# All data are saved as GeoTIFFs or ESRI Shapefiles in a specified directory.
# Data retrieved using the NLDI, Py3DEP, PyGeoHydro, and PyNHD packages.
from pathlib import Path

import py3dep as py3d
import pynhd as nhd
from pynhd import NLDI, NHDPlusHR, WaterData
import pygeohydro as gh
from pygeohydro import NWIS
import geopandas as gpd
import matplotlib.pyplot as plt
import warnings

# %%
# storage path
storage_path  = Path("/storage/dlhogan/summa_modeling_data/")
# basin name and outlet gauge id
name = "domain_EastRiver"
if name == "domain_EastRiver":
    gauge_id = "09112500"
elif name == "domain_TuolumneRiver":
    gauge_id = "11274790"
data_path = storage_path / name
# crs to convert to 
data_crs = 5070

# %%
def get_basin_geom(name=name, gauge_id=gauge_id, data_path=data_path, save=True):
    """
    Retrieve the basin geometry for a given USGS gauge site using the NLDI (Network-Linked Data Index) service.

    Parameters:
    -----------
    name : str
        The name to use when saving the basin shapefile.
    gauge_id : str
        The USGS gauge site identifier.
    data_path : Path
        The directory path where the basin geometry will be saved (if save=True).
    save : bool, optional (default=True)
        If True, saves the retrieved basin geometry as an ESRI Shapefile in the specified data_path.

    Returns:
    --------
    tuple
        basin : GeoDataFrame
            A GeoDataFrame containing the basin geometry.
        basin_geometry : geometry-like
            The extracted basin geometry as a shapely Polygon or MultiPolygon.

    Notes:
    ------
    - The function retrieves the basin geometry using the `NLDI.get_basins` method.
    - If `save=True`, the geometry is saved as a shapefile named `{name}.shp` in the given directory.
    """

    nldi = NLDI()
    
    basin = nldi.get_basins(gauge_id)
    # save the geometry to a file
    if save:
        basin.to_file(data_path / "shapefiles" / "catchment" / f"{name}.shp", driver="ESRI Shapefile")

    basin_geometry = basin.geometry.iloc[0]
    return basin, basin_geometry

# %%
def get_topo(basin_geometry, data_path=data_path, save=True):
    """
    Retrieve topographic data for a given basin, including elevation, slope, and aspect.

    Parameters:
    -----------
    basin_geometry : geometry-like
        The basin geometry (e.g., a shapely polygon) defining the area of interest.
    data_path : Path
        The directory path where the topographic data will be saved (if save=True).
    save : bool, optional (default=True)
        If True, saves the retrieved topographic data as GeoTIFF files in the specified data_path.

    Returns:
    --------
    topo : dict
        A dictionary containing topographic datasets for the specified basin, including:
        - 'elevation': Digital Elevation Model (DEM).
        - 'slope_degrees': Slope in degrees.
        - 'aspect_degrees': Aspect in degrees.

    Notes:
    ------
    - The function retrieves topographic data using the `py3d.get_map` method.
    - The data is retrieved at a 30-meter resolution.
    - The source CRS is EPSG:4326, and the output CRS is EPSG:5070.
    - If `save=True`, the data is saved as GeoTIFF files in the given directory.
    """
    topo = py3d.get_map(["DEM", "Slope Degrees", "Aspect Degrees"], basin_geometry, 30, geo_crs=4326, crs=5070)

    if save:
        # save each individual layer to a geotiff 
        topo["elevation"].rio.to_raster(data_path / "parameters" / "dem" / "dem.tif")
        topo["slope_degrees"].rio.to_raster(data_path / "parameters" / "dem" / "slope.tif")
        topo["aspect_degrees"].rio.to_raster(data_path / "parameters" / "dem" / "aspect.tif")
    return topo

# %%
def get_land_cover(basin, data_path=data_path, data_crs=data_crs, save=True):
    """
    Retrieve land cover data for a given basin using the NLCD (National Land Cover Database).

    Parameters:
    -----------
    basin : geometry-like
        The basin geometry (e.g., a shapely polygon) defining the area of interest.
    data_path : Path
        The directory path where the land cover data will be saved (if save=True).
    data_crs : str or CRS
        The coordinate reference system (CRS) of the input basin geometry.
    save : bool, optional (default=True)
        If True, saves the retrieved land cover data as GeoTIFF files in the specified data_path.

    Returns:
    --------
    land_cover : dict
        A dictionary containing land cover datasets for the specified basin, including:
        - 'cover_2021': Land cover classification for 2021.
        - 'descriptor_2021': Urban descriptor classification for 2021.
        - 'canopy_2021': Tree canopy cover for 2021.

    Notes:
    ------
    - The function retrieves land cover data using the `gh.nlcd_bygeom` method.
    - If `save=True`, the data is saved as GeoTIFF files in the given directory.
    - The function extracts data for the year 2021.
    """
    # download land cover data from NLCD
    land_cover = gh.nlcd_bygeom(basin, crs=data_crs, years={'cover':2021, 'descriptor':2021, 'canopy':2021})
    if save:
        # save the land cover data to a geotiff
        land_cover[list(land_cover.keys())[0]]['cover_2021'].rio.to_raster(data_path / "parameters" / "landclass" / "nlcd_land_cover_2021.tif")
        land_cover[list(land_cover.keys())[0]]['descriptor_2021'].rio.to_raster(data_path / "parameters" / "landclass" / "nlcd_urban_descriptor_2021.tif")
        land_cover[list(land_cover.keys())[0]]['canopy_2021'].rio.to_raster(data_path / "parameters" / "landclass" / "nlcd_canopy_2021.tif")
    return land_cover

# %%
def get_streamflow(gauge_id=gauge_id, data_path=data_path, freq='dv', dates=("2021-10-01","2023-09-30"), save=True):
    """
    Retrieve streamflow observations for a given USGS gauge site from the NWIS (National Water Information System).

    Parameters:
    -----------
    gauge_id : str
        The USGS gauge site identifier.
    data_path : Path
        The directory path where the streamflow data will be saved (if save=True).
    freq : str, optional (default='dv')
        The frequency of the streamflow data. Options include:
        - 'dv' for daily values
        - 'iv' for instantaneous values
    dates : tuple of str, optional (default=("2021-10-01", "2023-09-30"))
        The date range for retrieving streamflow data, formatted as (start_date, end_date).
    save : bool, optional (default=True)
        If True, saves the retrieved streamflow data as a CSV file in the specified data_path.

    Returns:
    --------
    q_obs : DataFrame
        A pandas DataFrame containing the retrieved streamflow observations.

    Notes:
    ------
    - The function retrieves streamflow data using the NWIS API.
    - If `save=True`, the data is saved as "streamflow_obs.csv" in the given directory.
    """
    # download streamflow observations
    q_obs = NWIS().get_streamflow(gauge_id, freq=freq, dates=dates, mmd=True)
    if save:
        # save the streamflow observations to a csv
        q_obs.to_csv(data_path / "benchmark" / "streamflow_obs.csv", index=True)
    return q_obs

# %%
def get_flowlines(gauge_id=gauge_id, data_path=data_path, save=True):    
    """
    Retrieve upstream flowlines for a given USGS gauge site, including both main 
    stem and tributary flowlines, using the NLDI (Network-Linked Data Index) service.

    Parameters:
    -----------
    gauge_id : str
        The USGS gauge site identifier.
    data_path : Path
        The directory path where the flowline shapefiles will be saved (if save=True).
    save : bool, optional (default=True)
        If True, saves the retrieved flowlines as ESRI Shapefiles in the specified data_path.

    Returns:
    --------
    tuple of GeoDataFrames
        flw_main : GeoDataFrame
            The main stem upstream flowlines.
        flw_trib : GeoDataFrame
            The upstream tributary flowlines.

    Notes:
    ------
    - The function retrieves flowlines up to 1000 km upstream.
    - Suppresses warnings related to long column names in ESRI Shapefiles.
    """
    # get the flowlines for the main stem and tributaries
    flw_main = NLDI().navigate_byid(
        fsource="nwissite",
        fid=f"USGS-{gauge_id}",
        navigation="upstreamMain",
        source="flowlines",
        distance=1000,
    )

    flw_trib = NLDI().navigate_byid(
        fsource="nwissite",
        fid=f"USGS-{gauge_id}",
        navigation="upstreamTributaries",
        source="flowlines",
        distance=1000,
    )
    if save:
        # ignore UserWarning UserWarning: Column names longer than 10 characters will be truncated when saved to ESRI Shapefile.
        warnings.filterwarnings("ignore", message="Column names longer than 10 characters will be truncated when saved to ESRI Shapefile.")
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        # save the flowlines to a file
        flw_main.to_file(data_path / "shapefiles" / "river_network" "flowlines_main.shp", driver="ESRI Shapefile")
        flw_trib.to_file(data_path / "shapefiles" / "river_network" "flowlines_trib.shp", driver="ESRI Shapefile")
    return flw_main, flw_trib

# %%
def get_subcatchments(basin, data_path=data_path, save=True):
    """
    Retrieve and optionally save subcatchments that intersect with a given basin.

    Parameters:
    basin (GeoDataFrame): A GeoDataFrame containing the basin geometry.
    data_path (Path, optional): The file path where the subcatchments shapefile will be saved. Defaults to a predefined data_path.
    save (bool, optional): If True, the subcatchments will be saved to a shapefile. Defaults to True.

    Returns:
    GeoDataFrame: A GeoDataFrame containing the subcatchments that intersect with the basin.
    """
    # gather all subcatchments
    wd_cat = WaterData("catchmentsp")
    cat = wd_cat.bygeom(basin.geometry.iloc[0], predicate="intersects")
    subcatchments = cat.clip(basin.geometry, )

    if save:
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        # save subcatchments to a shapefile
        subcatchments.to_file(data_path / "shapefiles" / "catchment" / "subcatchments.shp", driver="ESRI Shapefile")
    return subcatchments

# %%


# %%
if __name__ == "__main__":
    save=False
    basin, basin_geometry = get_basin_geom(save=save)
    topo = get_topo(basin_geometry, save=save)
    land_cover = get_land_cover(basin, save=save)
    q_obs = get_streamflow(save=True)
    flowlines = get_flowlines(save=save)
    subcatchments = get_subcatchments(basin,save=save)


