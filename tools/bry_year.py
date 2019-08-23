
import netCDF4

filename = "/lustre/storeB/project/fou/hi/arktis2030/metroms_run/barents-2.5km/BRY_2019.nc"

ds = netCDF4.Dataset(filename, mode="r+")

fill = dict()

for name, var in ds.variables.items():
    if "days" in var.dimensions and name != "Time":
        fill[name] = var[0,:]

for n in range(365):
    print("Handling day {}".format(n))
    
    ds.variables["Time"][n] = n

    for name, var in fill.items():
        ds.variables[name][n,:] = var[:]

ds.close()
