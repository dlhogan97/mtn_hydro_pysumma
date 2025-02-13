import cdsapi
import numpy as np

def download_ERA5_surface_pressure(bounding_box, years):
    """
    Downloads the ERA5 surface pressure data for the specified bounding box.
    
    Parameters:
        bounding_box (list): A list containing the coordinates of the bounding box in the format [north, west, south, east].
        years (list): A list of years for which to download the data.
    """
    dataset = "reanalysis-era5-single-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": ["surface_pressure"],
        "year": years,
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "day": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31"
        ],
        "time": [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ],
        "data_format": "netcdf",
        "download_format": "zip",
        "area": bounding_box
    }
    return dataset, request

tuolumne_bounding_box = [39, -120, 37, -119]
tuolumne_target = "/storage/dlhogan/summa_modeling_data/domain_TuolumneRiver/forcing/1_raw_data/ERA5_surface_pressure_2011_2024.nc"
east_river_bounding_box = [40, -107, 38, -105]
east_target = "/storage/dlhogan/summa_modeling_data/domain_EastRiver/forcing/1_raw_data/ERA5_surface_pressure_2011_2024.nc"
years = [str(v) for v in np.arange(2011, 2024).tolist()]

client = cdsapi.Client()

for i in range(2):
    if i == 0:
        bounding_box = tuolumne_bounding_box
        target = tuolumne_target
    else:
        bounding_box = east_river_bounding_box
        target = east_target
    dataset, request = download_ERA5_surface_pressure(tuolumne_bounding_box, years)
    client.retrieve(dataset, request, target).download()

# build a test suite for this function
def test_download_ERA5_surface_pressure():
    """
    Tests the download_ERA5_surface_pressure function.
    """
    bounding_box = [40, -107, 38, -105]
    years = ["2011", "2012"]
    dataset, request = download_ERA5_surface_pressure(bounding_box, years)
    
    assert dataset == "reanalysis-era5-single-levels"
    assert request["product_type"] == ["reanalysis"]
    assert request["variable"] == ["surface_pressure"]
    assert request["year"] == years
    return

