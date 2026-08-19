# Copernicus climate data (ERA5-Land)

The file in this folder (`data_stream-moda.nc`) is committed directly
rather than left as a download stub, since it's small enough to include.

- Source: Copernicus Climate Data Store (CDS), ERA5-Land monthly averaged
  reanalysis (`reanalysis-era5-land-monthly-means`):
  https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means
  
- Variables: 2m temperature (`2m_temperature`) and total precipitation
  (`total_precipitation`), monthly, 2000-2016, all 12 months.
  
- Area: Europe bounding box (North 71.5, West -12.5, South 34.0, East 35.0).

- Format: NetCDF.

- How it was obtained:*downloaded via the CDS API using the R `ecmwfr`
  package. The exact request (commented, since it needs your own CDS
  account credentials) is in `4_Fixed_effects_calculation.R`, in the
  `##CLIMATE` section - uncomment it and fill in your own `user`/`key` to
  reproduce the download from scratch.
  
- Units: the script converts temperature from Kelvin to Celsius and
  precipitation from meters to millimeters before extraction.
