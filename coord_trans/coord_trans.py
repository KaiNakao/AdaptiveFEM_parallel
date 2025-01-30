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

# vector (dE, dN, dU) to (dx, dy, dz)
def enu_to_xyz(lat, lon, de, dn, du, lat_c, lon_c):
    lat_rad = lat / 180 * np.pi
    lon_rad = lon / 180 * np.pi
    lat_c_rad = lat_c / 180 * np.pi
    lon_c_rad = lon_c / 180 * np.pi
    mat_c = np.array([[-np.sin(lon_c_rad), np.cos(lon_c_rad), 0],
                      [-np.sin(lat_c_rad) * np.cos(lon_c_rad), -np.sin(lat_c_rad)
                       * np.sin(lon_c_rad), np.cos(lat_c_rad)],
                      [np.cos(lat_c_rad) * np.cos(lon_c_rad), np.cos(lat_c_rad) * np.sin(lon_c_rad), np.sin(lat_c_rad)]])
    mat = np.array([[-np.sin(lon_rad), np.cos(lon_rad), 0],
                    [-np.sin(lat_rad) * np.cos(lon_rad), -np.sin(lat_rad)
                     * np.sin(lon_rad), np.cos(lat_rad)],
                    [np.cos(lat_rad) * np.cos(lon_rad), np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)]])
    return np.dot(np.linalg.inv(mat_c).T, np.dot(np.linalg.inv(mat), np.array([de, dn, du])))

# octree level
# nk = 7
nk = 1

# target moment tensor (miyagi)
target_lat = 38.4
target_lon = 141.2
target_depth = 25000.0
target_mxx = 0.0
target_myy = 0.0
target_mzz = 0.0
target_mxy = 0.0
target_myz = 0.0
target_mzx = 1e19
target_mvec = np.array([target_mxx, target_myy, target_mzz, target_mxy, target_myz, target_mzx])

# target area
# min_lon = 140.5
# max_lon = 143.5
min_lon = target_lon - 1.5
max_lon = target_lon + 1.5

# min_lat = 41.7
# max_lat = 43.7
min_lat = target_lat - 1.0
max_lat = target_lat + 1.0

# grid size
ds = 1250.


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
    # jivsm_region[hlabel] = jivsm_region[elvlabel] + geoid
    jivsm_region[hlabel] = jivsm_region[elvlabel]

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
        # h = elv + geoid
        h = elv
        layer[j] = lonlat_to_local(lat, lon, h, lat_c, lon_c, h_c)
    layers.append(layer)

# coordinate transformation of the observation data
# GNSS
column_names = ["lon", "lat", "elv", "dE", "dN", "dU", "id"]
df_gnss = pd.read_csv("../coord_trans/original/GNSS_org.dat", delim_whitespace=True, skiprows=1, names=column_names)

# exclude invalid observation
min_lon_ = min_lon + (max_lon - min_lon) * 0.25
max_lon_ = max_lon - (max_lon - min_lon) * 0.25
min_lat_ = min_lat + (max_lat - min_lat) * 0.25
max_lat_ = max_lat - (max_lat - min_lat) * 0.25
df_gnss = df_gnss[df_gnss["lon"] > min_lon_]
df_gnss = df_gnss[df_gnss["lon"] < max_lon_]
df_gnss = df_gnss[df_gnss["lat"] > min_lat_]
df_gnss = df_gnss[df_gnss["lat"] < max_lat_]
df_gnss.reset_index(inplace=True, drop=True)

# calculate geoid
lat = df_gnss["lat"]
lon = df_gnss["lon"]
with open("../gsigeo2011_ver2_1_asc/program/input.dat", "w") as f:
    for i in range(len(lat)):
        j = (i + 1) % 10000
        num = f'{j:4}'
        name = f'{j:18}'
        lat_str = f'{to_dms(lat[i]):15.04f}'
        lon_str = f'{to_dms(lon[i]):15.04f}'
        f.write(num + name + lat_str + lon_str + "\n")
subprocess.run(["./gsigeome_asc", "input.dat", "output.dat", "gsigeo2011_ver2_1.asc"],
                cwd="../gsigeo2011_ver2_1_asc/program")
geoid = np.loadtxt(
    "../gsigeo2011_ver2_1_asc/program/output.dat")[:, 4]
df_gnss["geoid"] = geoid
df_gnss["h"] = df_gnss["elv"] + df_gnss["geoid"]
df_gnss.to_csv("data/GNSS_reduced.dat", index=False, sep=" ")

df_gnss_out = pd.DataFrame()

# coordinate transformation for location
xvec = []
yvec = []
for i in df_gnss.index:
    lat = df_gnss.loc[i, "lat"]
    lon = df_gnss.loc[i, "lon"]
    h = df_gnss.loc[i, "elv"] + df_gnss.loc[i, "geoid"]
    xyz = lonlat_to_local(lat, lon, h, lat_c, lon_c, h_c)
    for j in range(3):
        xvec.append(xyz[0])
        yvec.append(xyz[1])

df_gnss_out["x"] = np.array(xvec)
df_gnss_out["y"] = np.array(yvec)

# observation direction
exvec = []
eyvec = []
ezvec = []
for i in df_gnss.index:
    d = [1, 0, 0]
    exvec.append(d[0])
    eyvec.append(d[1])
    ezvec.append(d[2])
    d = [0, 1, 0]
    exvec.append(d[0])
    eyvec.append(d[1])
    ezvec.append(d[2])
    d = [0, 0, 1]
    exvec.append(d[0])
    eyvec.append(d[1])
    ezvec.append(d[2])

df_gnss_out["ex"] = exvec
df_gnss_out["ey"] = eyvec
df_gnss_out["ez"] = ezvec

# line of sight displacement
dlosvec = []
for i in df_gnss.index:
    lat = df_gnss.loc[i, "lat"]
    lon = df_gnss.loc[i, "lon"]
    de = df_gnss.loc[i, "dE"]
    dn = df_gnss.loc[i, "dN"]
    du = df_gnss.loc[i, "dU"]
    disp = enu_to_xyz(lat, lon, de, dn, du, lat_c, lon_c)
    dlosvec.append(disp[0])
    dlosvec.append(disp[1])
    dlosvec.append(disp[2])
df_gnss_out["dlos"] = dlosvec

# observation error
sigmavec = []
for i in df_gnss.index:
    sigma = np.array([0.2, 0.3, 0.9])
    # increase error for incomplete data
    if df_gnss.loc[i, "id"] in ["950132", "950141"]:
        sigma *= 10
    sigmavec.append(sigma[0])
    sigmavec.append(sigma[1])
    sigmavec.append(sigma[2])
df_gnss_out["sigma"] = sigmavec

# type
df_gnss_out["type"] = ["GNSS"] * len(df_gnss_out)
df_gnss_out.set_axis(["x", "y", "ex", "ey", "ez", "dlos", "sigma", "type"], axis='columns')
df_gnss_out.to_csv("data/observation_gnss.dat", index=False, sep=" ")

# df_obs = pd.concat([df_sar_out, df_gnss_out])
df_obs = pd.concat([df_gnss_out])
df_obs.to_csv("data/observation.dat", index=False, sep=" ")

# coordinate of centroid
lat = [target_lat]
lon = [target_lon]
with open("../gsigeo2011_ver2_1_asc/program/input.dat", "w") as f:
    for i in range(len(lat)):
        j = (i + 1) // 10000
        num = f'{j:4}'
        name = f'{j:18}'
        lat_str = f'{to_dms(lat[i]):15.04f}'
        lon_str = f'{to_dms(lon[i]):15.04f}'
        f.write(num + name + lat_str + lon_str + "\n")
subprocess.run(["./gsigeome_asc", "input.dat", "output.dat", "gsigeo2011_ver2_1.asc"],
                cwd="../gsigeo2011_ver2_1_asc/program")
geoid = np.loadtxt(
    "../gsigeo2011_ver2_1_asc/program/output.dat")[4]
target_geoid = geoid
target_h = -target_depth + target_geoid
target_xyz = lonlat_to_local(target_lat, target_lon, target_h, lat_c, lon_c, h_c)

# find max/min of x and y
xmin = int(min(min(layer[:,0]) for layer in layers))
# xmax = int(max(max(layer[:,0]) for layer in layers))
xmax = xmin + 240000
ymin = int(min(min(layer[:,1]) for layer in layers))
# ymax = int(max(max(layer[:,1]) for layer in layers))
ymax = ymin + 200000

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

# deletion thin layers
for i in range(layers_new[14].shape[0]):
    layers_new[1][i, 2] = layers_new[0][i, 2]
    layers_new[2][i, 2] = layers_new[0][i, 2]
    layers_new[3][i, 2] = layers_new[0][i, 2]
    layers_new[4][i, 2] = layers_new[0][i, 2]
    layers_new[5][i, 2] = layers_new[0][i, 2]
    layers_new[6][i, 2] = layers_new[0][i, 2]
    layers_new[8][i, 2] = layers_new[7][i, 2]
    layers_new[9][i, 2] = layers_new[7][i, 2]
    layers_new[10][i, 2] = layers_new[7][i, 2]
    layers_new[11][i, 2] = layers_new[7][i, 2]
    layers_new[12][i, 2] = layers_new[7][i, 2]
    layers_new[13][i, 2] = layers_new[7][i, 2]
for i in range(layers[14].shape[0]):
    layers[1][i, :] = layers[0][i, :]
    layers[2][i, :] = layers[0][i, :]
    layers[3][i, :] = layers[0][i, :]
    layers[4][i, :] = layers[0][i, :]
    layers[5][i, :] = layers[0][i, :]
    layers[6][i, :] = layers[0][i, :]
    layers[8][i, :] = layers[7][i, :]
    layers[9][i, :] = layers[7][i, :]
    layers[10][i, :] = layers[7][i, :]
    layers[11][i, :] = layers[7][i, :]
    layers[12][i, :] = layers[7][i, :]
    layers[13][i, :] = layers[7][i, :]

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
    f.write(str(nk) + "\n")
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
        f.write(str(1000) + "\n")
        for _ in range(2):
            f.write(str(0.) + "\n")
    
with open("data/para_setting.dat", "w") as f:
    f.write("num of MPI processes\n")
    f.write(str(1)+"\n")
    f.write("num of OpenMP threads\n")
    f.write(str(1)+"\n")

with open("data/log_setting.dat", "w") as f:
    f.write("output log for inner CGs (1:on)\n")
    f.write(str(0) + "\n")
    f.write("output summary for inner CG per outer iteration (1:on)\n")
    f.write(str(0) + "\n")

# pointload
with open("data/pointload.dat", "w") as f:
    f.write("number of loads\n")
    f.write(str(df_obs.shape[0]) + "\n")
    f.write("x y ex ey ez\n")
    xvec = df_obs["x"].values
    yvec = df_obs["y"].values
    exvec = df_obs["ex"].values
    eyvec = df_obs["ey"].values
    ezvec = df_obs["ez"].values
    for iobs in range(df_obs.shape[0]):
        x = xvec[iobs] - xmin
        y = yvec[iobs] - ymin
        ex = exvec[iobs]
        ey = eyvec[iobs]
        ez = ezvec[iobs]
        f.write("%f %f %f %f %f\n" % (x, y, ex, ey, ez))

# target moment tensor
with open("data/target_centroid.dat", "w") as f:
    f.write("centroid coordinate\n")
    #f.write(str(target_xyz[0] - xmin) + "\n")
    #f.write(str(target_xyz[1] - ymin) + "\n")
    #f.write(str(target_xyz[2] - zmin) + "\n")
    f.write(str(62000) + "\n")
    f.write(str(102000) + "\n")
    f.write(str(162000) + "\n")
    f.write("moment tensor\n")
    for i in range(6):
        f.write(str(target_mvec[i]) + "\n")
