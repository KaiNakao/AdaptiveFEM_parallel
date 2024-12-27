import numpy as np
import pandas as pd
import subprocess
from scipy.spatial import KDTree

def to_dms(value):
    d = int(value)
    value -= d
    m = int(value * 60)
    value -= m / 60
    s = (value * 3600)
    if 60 - s < 0.00001:
        m += 1
        s = 0
    if 60 - m < 0.00001:
        d += 1
        m = 0
    ret = d * 10000 + m * 100 + s
    return ret

def lonlat_to_ecef(lat, lon, h):
    lat_rad = lat / 180. * np.pi
    lon_rad = lon / 180. * np.pi
    a = 6378137
    e = 0.081819191042815791
    n = a / np.sqrt(1 - e**2 * np.sin(lat_rad)**2)
    x = (n + h) * np.cos(lat_rad) * np.cos(lon_rad)
    y = (n + h) * np.cos(lat_rad) * np.sin(lon_rad)
    z = (n * (1 - e**2) + h) * np.sin(lat_rad)
    return np.array([x, y, z])

def lonlat_to_local(lat, lon, h, lat_c, lon_c, h_c):
    ecef_c = lonlat_to_ecef(lat_c, lon_c, h_c)
    ecef = lonlat_to_ecef(lat, lon, h)
    lat_c_rad = lat_c / 180 * np.pi
    lon_c_rad = lon_c / 180 * np.pi
    mat = np.array([[-np.sin(lon_c_rad), np.cos(lon_c_rad), 0],
                    [-np.sin(lat_c_rad) * np.cos(lon_c_rad), -np.sin(lat_c_rad)
                     * np.sin(lon_c_rad), np.cos(lat_c_rad)],
                    [np.cos(lat_c_rad) * np.cos(lon_c_rad), np.cos(lat_c_rad) * np.sin(lon_c_rad), np.sin(lat_c_rad)]])
    return np.dot(mat, ecef - ecef_c)

# target area
min_lon = 141.3
max_lon = 142.6
min_lat = 42.2
max_lat = 43.2
# min_lon = 140.0
# max_lon = 144.0
# min_lat = 40.0
# max_lat = 44.0

# grid size
# ds = 1000.
ds = 2500.
# ds = 10000.

# number of layers for JIVSM
nlayer = 23

# number of points considered for near-neighbor interpolation
nneighbor = 4

# Extract JIVSM in target area
jivsm = np.loadtxt("../jivsm/Ejapan_path20111110.dat")
header = ["lon", "lat"] + ["elv" + str(i+1).zfill(2) for i in range(nlayer)]
jivsm = pd.DataFrame(jivsm, columns=header)
jivsm_region = jivsm[(min_lon <= jivsm["lon"]) & (jivsm["lon"] <= max_lon) & (min_lat <= jivsm["lat"]) & (jivsm["lat"] <= max_lat)]
jivsm_region.reset_index(inplace=True, drop=True)
jivsm_region.to_csv("../jivsm/Ejapan_path20111110_area.dat", index=False, sep=" ")

# calculate geoid for grid points of JIVSM
lat = jivsm_region.loc[:, "lat"]
lon = jivsm_region.loc[:, "lon"]
with open("../gsigeo2011_ver2_1_asc/program/input.dat", "w") as f:
    for i in range(jivsm_region.shape[0]):
        j = (i + 1) // 10000
        num = f'{j:4}'
        name = f'{j:18}'
        lat_str = f'{to_dms(lat[i]):15.04f}'
        lon_str = f'{to_dms(lon[i]):15.04f}'
        f.write(num + name + lat_str + lon_str + "\n")
subprocess.run(["./gsigeome_asc", "input.dat", "output.dat", "gsigeo2011_ver2_1.asc"],
                cwd="../gsigeo2011_ver2_1_asc/program")
geoid = np.loadtxt(
    "../gsigeo2011_ver2_1_asc/program/output.dat")[:, 4]
jivsm_region["geoid"] = geoid

# calculate ellipsoidal height (elevation + geoid)
for i in range(nlayer):
    hlabel = "h" + str(i+1).zfill(2)
    elvlabel = "elv" + str(i+1).zfill(2)
    jivsm_region[hlabel] = jivsm_region[elvlabel] + geoid

# define origin of cartesian coordinates
max_h = max(jivsm_region["h01"])
min_h = min(jivsm_region["h23"])
lat_c = (min_lat + max_lat) / 2.
lon_c = (min_lon + max_lon) / 2.
h_c = (min_h + max_h) / 2.

# coordinate transformation for each layer
layers = []
for i in range(nlayer):
    layer = np.zeros((jivsm_region.shape[0], 3))
    for j in range(jivsm_region.shape[0]):
        lat = jivsm_region.loc[j, "lat"]
        lon = jivsm_region.loc[j, "lon"]
        geoid = jivsm_region.loc[j, "geoid"]
        elv = jivsm_region.loc[j, "elv"+str(i+1).zfill(2)]
        h = elv + geoid
        layer[j] = lonlat_to_local(lat, lon, h, lat_c, lon_c, h_c)
    layers.append(layer)

# find max/min of x and y
xmin = int(min(min(layer[:,0]) for layer in layers))
xmax = int(max(max(layer[:,0]) for layer in layers))
ymin = int(min(min(layer[:,1]) for layer in layers))
ymax = int(max(max(layer[:,1]) for layer in layers))

# number of intervals for new grid in cartesian coordinates
nx = int((xmax - xmin) / ds)
ny = int((ymax - ymin) / ds)

# grid point
points_new = []
for iy in range(ny + 1):
    y = ymin + iy * ds
    for ix in range(nx + 1):
        x = xmin + ix * ds
        points_new.append([x, y])
points_new = np.array(points_new)

# interpolate zvalue for new grid point from old grid
layers_new = []
for layer in layers:
    tree = KDTree(layer[:,:2])
    dd, ii = tree.query(points_new, k=4)
    rr = dd**(-1)
    zz = layer[:, 2]
    layer_new = np.zeros((points_new.shape[0], 3))
    for i in range(layer_new.shape[0]):
        layer_new[i, 0] = points_new[i, 0]
        layer_new[i, 1] = points_new[i, 1]
        layer_new[i, 2] = np.dot(rr[i], zz[ii[i]]) / sum(rr[i])
    layers_new.append(layer_new)

# minimum of x, y, z is zero in meshing
zmax = min(layers_new[0][:,2])
zmin = min(layers_new[-1][:,2])
zmin -= 0.2 * (zmax - zmin)
for i in range(nlayer):
    layers_new[i][:, 0] -= xmin
    layers_new[i][:, 1] -= ymin
    layers_new[i][:, 2] -= zmin

# output
for i in range(nlayer):
    filename = "data/sur" + str(i+1).zfill(4) + ".dat"
    with open(filename, "w") as f:
        f.write("nx, ny\n")
        f.write(str(nx) + " " + str(ny) + "\n")
        f.write("DEM data\n")
        for j in range(layer_new.shape[0]):
            f.write(str(layers_new[i][j, 2]) + "\n")

with open("data/modeling_setting.dat", "w") as f:
    f.write("nx, ny\n")
    f.write(str(nx) + " " + str(ny) + "\n")
    f.write("ds\n")
    f.write(str(ds) + "\n")
    f.write("num of layer\n")
    f.write(str(nlayer) + "\n")
    f.write("nk 2**nk: nx and ny = 2**nk * alpha\n")
    f.write(str(4) + "\n")
    f.write("freq\n")
    f.write(str(0.25) + "\n")
    f.write("now\n")
    f.write(str(5) + "\n")

material = np.loadtxt("../jivsm/material.dat", skiprows=1)
with open("data/material.dat", "w") as f:
    for i in range(nlayer):
        f.write("matrial " + str(i + 1) + "\n")
        f.write(str(material[i, 0] * 1000) + "\n")
        f.write(str(material[i, 1] * 1000) + "\n")
        f.write(str(material[i, 2] * 1000) + "\n")
        for _ in range(4):
            f.write(str(0.) + "\n")
        f.write(str(1.) + "\n")
        for _ in range(2):
            f.write(str(0.) + "\n")
    
with open("data/para_setting.dat", "w") as f:
    f.write("num of MPI processes\n")
    f.write(str(1)+"\n")
    f.write("num of OpenMP threads\n")
    f.write(str(1)+"\n")

with open("data/log_setting.dat", "w") as f:
    f.write("output log for inner CGs (1:on)\n")
    f.write(str(1) + "\n")
    f.write("output summary for inner CG per outer iteration (1:on)\n")
    f.write(str(1) + "\n")

# fault_x = 16500.
# fault_y = 20500.
# fault_z = 17500.
# strike = 90.
# dip = 90.
# rake = 0.
# moment = 1.0 * 10**19
# with open("data/faultpara.dat", "w") as f:
#     f.write("num of point source\n")
#     f.write("1\n")
#     f.write("1 th point source\n")
#     f.write(str(fault_x) + "\n")
#     f.write(str(fault_y) + "\n")
#     f.write(str(fault_z) + "\n")
#     f.write(str(strike) + "\n")
#     f.write(str(dip) + "\n")
#     f.write(str(rake) + "\n")
#     f.write(str(moment) + "\n")

# pointload (tentative)
with open("data/pointload.dat", "w") as f:
    x = (xmin + xmax) / 2.
    y = (ymin + ymax) / 2.
    ex = -0.5900865906800314 
    ey = 0.10467736655876937 
    ez = -0.8005251179256889
    f.write("x y ex ey ez\n")
    f.write("%f %f %f %f %f\n" % (x, y, ex, ey, ez))
