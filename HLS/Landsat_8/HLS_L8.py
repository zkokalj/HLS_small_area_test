import os
import rasterio
from rasterio import warp
from xml.dom import minidom
import numpy as np
from scipy import ndimage
from scipy.interpolate import lagrange
import numpy.polynomial.polynomial as nppol
from numpy.polynomial.polynomial import Polynomial
import time
import shutil
from openpyxl import load_workbook
from math import floor

"""
Process HLS algorithm on Landsat-8 images. One directory is mandatory : result of LaSRC algorithm. 
One xlsx file is required : Landsat Collection 1 to Collection 2 PR Mean Offsets (https://www.usgs.gov/media/files/landsat-collection-1-collection-2-pr-mean-offsets)
"""

#Path of initial image
input_initial = "C1L1"

#Path of LaSRC result
input_LaSRC_result = "LaSRC_result"

#Path of result images
output_image = "Result"

#path of Landsat Collection 1 to Collection 2 PR Mean Offsets file
path_xlsx = "C1_C2_WRS-2_Path_Row_Mean_Offsets.xlsx"



#Polynomial coefficients to calculate normalized sun zenithal angle
table_theta_s_out = [31, -0.127, 0.0119, 2.39*10**(-5), -9.46*10**(-7), -1.94*10**(-9), 6.13*10**(-11)]

#Values of BRDF-coefficients
f_iso = [0.0774, 0.774, 0.1306, 0.1690, 0.3093, 0.3430, 0.2658]
f_geo = [0.0079, 0.079, 0.0178, 0.0227, 0.0330, 0.0453, 0.0387]
f_vol = [0.0372, 0.0372, 0.0580, 0.0574, 0.1535, 0.1154, 0.0639]

#Polynomes for Fmask processing
x = np.array([0, 1, 2, 3, 4, 255])
y = np.array([0, 32, 8, 16, 2, 255])
convert_Fmask_Sentinel = lagrange(x, y)
convert_Fmask_Sentinel=list(Polynomial(convert_Fmask_Sentinel).coef)
convert_Fmask_Sentinel.reverse()

x1 = np.array([0, 32, 8, 16, 2, 255])
y1 = np.array([0, 0, 1, 0, 1, 0])
cloud_shadow_Sentinel = lagrange(x1, y1)
cloud_shadow_Sentinel=list(Polynomial(cloud_shadow_Sentinel).coef)
cloud_shadow_Sentinel.reverse()

x2 = np.array([0, 2, 4, 8, 16, 32, 255])
y2 = np.array([0, 1, 2, 3, 4, 5, 6])
small_values = lagrange(x2, y2)
small_values=list(Polynomial(small_values).coef)
small_values.reverse()

#Band in Landsat-8 images
list_band = [["_Fmask4.tif", 1], ["_sr_band1.tif", 1], ["_sr_band2.tif", 1], ["_sr_band3.tif", 1], ["_sr_band4.tif", 1], ["_sr_band5.tif", 1],  ["_sr_band6.tif", 1], ["_sr_band7.tif", 1], ["_toa_band9.tif", 1], ["_B10.tif", 2], ["_B11.tif", 2]]

#Spectral bands processed by BRDF normalization
list_band_brdf = ["_sr_band1.tif", "_sr_band2.tif", "_sr_band3.tif", "_sr_band4.tif", "_sr_band5.tif",  "_sr_band6.tif",  "_sr_band7.tif"]

#Spectral bands not processed by BRDF normalization
list_band_without_brdf = ["_toa_band9.tif", "_B10.tif", "_B11.tif"]

#Thermal bands
list_thermal_band = [["_B10.tif", "10"], ["_B11.tif", "11"]]

def read_xlsx():
    """
    Read xlsx file which contains translation values for coregistration
    """
    wb = load_workbook(filename=path_xlsx, read_only=True)
    ws = wb['Mean Offsets']
    return ws

def find_translation(file, ws):
    """
    Find translation for coregistration
    """
    WRS_path = int(file[10:13])
    WRS_row = int(file[13:16])
    for row in ws.rows:
        if row[0].value==WRS_path and row[1].value==WRS_row:
            return row[2].value, row[3].value

def open_tif(path_file):
    """
    Open tif file and return values in radian
    """
    dataset = rasterio.open(path_file)
    array = dataset.read(1)
    array[array==-32768] = 0
    return array / 100 * np.pi / 180

def find_centre(path_file, rows, cols):
    """
    Find latitude of center of the image with information in MTL.txt file and return solar angle for BRDF-Normalization.
    """
    file_ANG_txt = open(path_file)
    somme_latitude = 0
    for line in file_ANG_txt:
        list_line = line.split(" ")
        if len(list_line) >= 7 and (list_line[4] == "CORNER_UL_LAT_PRODUCT" or list_line[4] == "CORNER_UR_LAT_PRODUCT" or list_line[4] == "CORNER_LL_LAT_PRODUCT" or list_line[4] == "CORNER_LR_LAT_PRODUCT"):
            somme_latitude += float(list_line[6][:-1])
    latitude =  somme_latitude / 4

    #Calculate polynome to get theta_s_out (solar angle for BRDF-Normalization)
    sum = 0
    for i in range(len(table_theta_s_out)):
        sum += table_theta_s_out[i] * latitude**i
    return np.ones((rows, cols)) * sum * np.pi / 180


def K_geo(t_v, t_s, phi):
    """
    Calculate K_geo kernel.
    t_v : zenithal angle view
    t_s : zenithal sun view
    phi : view-sun relative azimutal angle
    """

    tan_tv = np.tan(t_v)
    tan_ts = np.tan(t_s)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_tv = np.cos(t_v)
    cos_ts = np.cos(t_s)
    sec_tv = 1 / cos_tv
    sec_ts = 1 / cos_ts
    sin_tv = np.sin(t_v)
    sin_ts = np.sin(t_s)


    D = np.sqrt( tan_tv**2 +  tan_ts**2 - 2 * tan_tv * tan_ts * cos_phi )

    cos_t = 2 * np.sqrt( D**2 + (tan_tv * tan_ts * sin_phi)**2) / (sec_tv + sec_ts)

    if np.max(cos_t) > 1:
        cos_t = np.ones(cos_t.shape)
        print("Alert : cos(t) is higher than 1. There will be problems in the results")

    t = np.arccos(cos_t)

    O = 1 / np.pi * (t - np.sin(t) * cos_t) * (sec_tv + sec_ts)

    cos_epsilon = cos_tv * cos_ts + sin_tv * sin_ts * cos_phi

    result = O - sec_tv - sec_ts + 0.5 * (1 + cos_epsilon) * sec_tv * sec_ts
    return result

def K_vol(t_v, t_s, phi):
    """
    Calculate K_vol kernel.
    t_v : zenithal angle view
    t_s : zenithal sun view
    phi : view-sun relative azimutal angle
    """

    cos_tv = np.cos(t_v)
    cos_ts = np.cos(t_s)
    sin_tv = np.sin(t_v)
    sin_ts = np.sin(t_s)
    cos_phi = np.cos(phi)

    cos_epsilon = cos_tv * cos_ts + sin_tv * sin_ts * cos_phi

    epsilon = np.arccos(cos_epsilon)

    result = ( (np.pi / 2 - epsilon) * cos_epsilon + np.sin(epsilon)) / (cos_tv + cos_ts) - np.pi / 4
    return result


def process_Fmask(path_file, file, directory_main, dst_transform):
    """
    Process Fmask : modify values, add adjacent cloud classes, transform values into a bit index.
    """
    path_Fmask = os.path.join(path_file, file+"_Fmask4.tif")
    if os.path.exists(path_Fmask):
        dataset = rasterio.open(path_Fmask)

        #Modify values to get values in HLS specification
        Fmask = nppol.polyval(dataset.read(1), convert_Fmask_Sentinel)
        Fmask = (Fmask+0.1).astype("i8")

        #Build an array with 1 for cloud and cloud shadow, 0 for others pixels
        cloud_shadow = nppol.polyval(Fmask, cloud_shadow_Sentinel)
        cloud_shadow = (cloud_shadow+0.1).astype("i8")

        #Dilatation
        structure = ndimage.generate_binary_structure(2, 2)
        dilatation = ndimage.binary_dilation(cloud_shadow, structure=structure, iterations=5).astype(cloud_shadow.dtype)

        #Modify value after dilatation adding class adjacent cloud / cloud shadow
        Fmask2 = adjacent_cloud(Fmask, cloud_shadow, dilatation)

        channel_small_values = nppol.polyval(Fmask2, small_values)
        #Calculate translation in pixel between initial Fmask and destination layer
        delta_i = floor((- 15 - translation_latitude) / 30)
        delta_j = floor((translation_longitude + - 15) / 30)

        #Build filter for convolution
        window = create_window(delta_i, delta_j)

        #Convolution
        resample = ndimage.convolve(channel_small_values, window)
        resample = (resample+0.1).astype(np.int32)

        #Transform result of bit index into a convolution
        vFmask_create_bit_index = np.vectorize(convert_bit_index)
        Fmask_bit_index = vFmask_create_bit_index(resample)

        save(os.path.join(directory_main, file+"_Fmask.tif"), Fmask_bit_index, dst_transform, np.uint16, 255)


def adjacent_cloud(Fmask, cloud_shadow, dilatation):
    """
    Return Fmask with all classes (adjacent cloud included).
    """
    cloud_shadow_Fmask = cloud_shadow * Fmask
    dilatation_Fmask = (dilatation - cloud_shadow) * 4
    other_Mask = - (dilatation - 1)
    remain_Fmask = other_Mask * Fmask
    return cloud_shadow_Fmask + dilatation_Fmask + remain_Fmask

def create_window(delta_i, delta_j):
    """
    Create window for convolution
    """
    #calculate size of the window
    if delta_i < 0:
        delta_i_max = - delta_i
    else:
        delta_i_max = delta_i + 1

    if delta_j < 0:
        delta_j_max = - delta_j
    else:
        delta_j_max = delta_j + 1

    max_delta = max(delta_i_max, delta_j_max)
    window_size = 2 * max_delta + 1
    window = np.zeros((window_size, window_size))

    #Complete window
    window[max_delta + delta_i][max_delta + delta_j] = 1
    window[max_delta + delta_i][max_delta + delta_j + 1] = 10
    window[max_delta + delta_i + 1][max_delta + delta_j] = 100
    window[max_delta + delta_i + 1][max_delta + delta_j + 1] = 1000
    window = window.astype(np.int16)
    return window

def convert_bit_index(x):
    """
    Transform x into a bit index
    """
    x_string = str(x)
    sum = 0
    if x == 6666:#if it is a no-data pixel
        return 255
    for i in range(1,6):
        if str(i) in x_string:
            sum += 2**i
    return sum

def coregistration(directory_initial, directory_LaSRC, directory_workspace, translation_latitude, translation_longitude, file):
    """
    Process coregistration on image and return rasterio. Affine object for final position of L30 images.
    """
    for name in list_band:
        if name[1] == 1:
            dataset = rasterio.open(os.path.join(directory_LaSRC, file+name[0]))
        else:
            dataset = rasterio.open(os.path.join(directory_initial, file+name[0]))
        channel = dataset.read(1)
        dst_transform = dataset.transform
        #rasterio.Affine object for coregistration
        dst_transform_co = rasterio.Affine(dst_transform[0], dst_transform[1], dst_transform[2] + translation_longitude, dst_transform[3], dst_transform[4], dst_transform[5] + translation_latitude)
        #rasterio.Affine object for final position of L30 images
        dst_transform_final = rasterio.Affine(dst_transform[0], dst_transform[1], dst_transform[2] - 15, dst_transform[3], dst_transform[4], dst_transform[5] + 15)

        nan_array = channel
        nan_array = nan_array.astype(np.float32)
        nan_array[channel == -9999] = np.nan

        if name[1] == 2:
            nan_array[channel == 0] = np.nan

        if name[0] != "_Fmask4.tif" :
            path_output = os.path.join(directory_workspace, file+name[0]+"_intermediaire.tif")
            nan_array = nan_array / 10000
        else:
            path_output = os.path.join(directory_workspace, file+name[0])
        #Save spectral band
        with rasterio.open(
                path_output,
                'w',
                driver='GTiff',
                width=nan_array.shape[1],
                height=nan_array.shape[0],
                count=1,
                dtype=np.float32,
                nodata=np.nan,
                transform=dst_transform_co,
                crs="EPSG:32633") as dst:
            dst.write(nan_array, indexes=1)
    return dst_transform_final

def resample(directory_workspace, file, directory_main, dst_transform):
    """
    Resample each spectral band of an image. A cubic convolution is used
    """
    for name in list_band_brdf:
        #Adjust image on Sentinel-2 image with cubic convolution
        resample_band(os.path.join(directory_workspace, file+name+"_intermediaire.tif"), os.path.join(directory_workspace, file+name), dst_transform, 1, np.float32)
    for name in list_band_without_brdf:
        if name == "_toa_band9.tif":
        #Adjust image on Sentinel-2 image with cubic convolution
            resample_band(os.path.join(directory_workspace, file+name+"_intermediaire.tif"), os.path.join(directory_main, file+name), dst_transform, 10000, np.int16)
        else:
            resample_band(os.path.join(directory_workspace, file+name+"_intermediaire.tif"), os.path.join(directory_workspace, file+name), dst_transform, 1, np.float32)


def resample_band(input, output, dst_transform, factor, dtype):
    """
    Resample spectral band.
    """
    dataset = rasterio.open(input)
    destination = np.zeros((dataset.height, dataset.width), np.float32)
    rasterio.warp.reproject(
            dataset.read(1),
            destination,
            src_transform=dataset.transform,
            src_crs="EPSG:32633",
            dst_transform=dst_transform,
            dst_crs="EPSG:32633",
            resampling=rasterio.warp.Resampling.cubic)


    #Save spectral band.
    with rasterio.open(
                output,
                'w',
                driver='GTiff',
                width=dataset.width,
                height=dataset.height,
                count=1,
                dtype=dtype,
                nodata=0,
                transform=dst_transform,
                crs="EPSG:32633") as dst:
        dst.write(destination * factor, indexes=1)



def BRDF(directory_workspace, num_band, file, K_geo_sensor, K_geo_norm, K_vol_sensor, K_vol_norm, directory_main):
    """
    Process BRDF-normalization on band num_band.
    """
    path_band = os.path.join(directory_workspace, file+list_band_brdf[num_band])
    if os.path.exists(path_band):
        dataset = rasterio.open(path_band)
        channel = dataset.read(1)
        f_iso_l = f_iso[num_band]
        f_geo_l = f_geo[num_band]
        f_vol_l = f_vol[num_band]

        #Calculate c-factor
        c = ( f_iso_l + f_geo_l * K_geo_norm + f_vol_l * K_vol_norm) / ( f_iso_l + f_geo_l * K_geo_sensor + f_vol_l * K_vol_sensor )

        #Calculate BRDF-normalized image
        BRDF_image = c * channel
        save(os.path.join(directory_main, file+list_band_brdf[num_band]), BRDF_image * 10000, dataset.transform, np.int16, 0)

def thermal_band(directory_workspace, directory_main, num_band, file):
    path_band = os.path.join(directory_workspace, file+list_thermal_band[num_band][0])
    if os.path.exists(path_band):
        dataset = rasterio.open(path_band)
        channel = dataset.read(1)
        K1 = float(find_data_txt(os.path.join(directory_initial, file+"_MTL.txt"), "K1_CONSTANT_BAND_"+list_thermal_band[num_band][1], "float"))
        K2 = float(find_data_txt(os.path.join(directory_initial, file+"_MTL.txt"), "K2_CONSTANT_BAND_"+list_thermal_band[num_band][1], "float"))
        channel = channel.astype(np.float64)
        channel[channel==0]=np.nan
        array = 10000 * 0.0003342 * channel + 0.1
        T_K = K2 / np.log(K1 / array + 1)
        T_C = T_K - 273.15
        save(os.path.join(directory_main, file+"_bt_band"+list_thermal_band[num_band][1]+".tif"), T_C*100, dataset.transform, np.int16, 0)

def save(path, array, dst_transform, type, nodata):
    """
    Save array in path like a tif with projection EPSG:32633.
    """
    with rasterio.open(
            path,
            'w',
            driver='GTiff',
            width=array.shape[1],
            height=array.shape[0],
            count=1,
            dtype=type,
            nodata=nodata,
            transform=dst_transform,
            crs="EPSG:32633") as dst:
        dst.write(array, indexes=1)

def create_xml_metadata(directory_LaSRC, file, directory_main, solar_azimuth_angle, solar_zenithal_angle, angle_view_zenith, angle_view_azimuth, theta_s_out):
    """
    Create xml file with metadata
    """

    #Open Fmask layer
    Fmask = rasterio.open(os.path.join(directory_main, file+"_Fmask.tif")).read(1)
    n, m = Fmask.shape
    count_Fmask_None = np.count_nonzero(Fmask == 255)
    count_Fmask_cloud = np.count_nonzero(Fmask == 2)
    count_Fmask_cloud_shadow = np.count_nonzero(Fmask == 4)

    #Create xml file
    newdoc = minidom.Document()
    root = newdoc.createElement('root')
    newdoc.appendChild(root)

    #Complete xml file
    LANDSAT_PRODUCT_ID = find_data_txt(os.path.join(directory_LaSRC, file+"_MTL.txt"), "LANDSAT_PRODUCT_ID", "string")
    add_node_xml(newdoc, root, 'LANDSAT_PRODUCT_ID', LANDSAT_PRODUCT_ID)
    SENSING_TIME = find_data_txt(os.path.join(directory_LaSRC, file+"_MTL.txt"), "SCENE_CENTER_TIME", "string")
    add_node_xml(newdoc, root, 'SENSING_TIME', SENSING_TIME)
    add_node_xml(newdoc, root, 'SPATIAL_COVERAGE', str(100 - count_Fmask_None / (n*m) * 100))
    add_node_xml(newdoc, root, 'CLOUD_COVERAGE', str((count_Fmask_cloud + count_Fmask_cloud_shadow) / (n*m - count_Fmask_None) * 100))
    add_node_xml(newdoc, root, 'SPATIAL_RESAMPLING_ALG', "Cubic convolution")
    ULX = find_data_txt(os.path.join(directory_LaSRC, file+"_MTL.txt"), "CORNER_UL_LAT_PRODUCT", "float")
    add_node_xml(newdoc, root, 'ULX', ULX)
    ULY = find_data_txt(os.path.join(directory_LaSRC, file+"_MTL.txt"), "CORNER_UL_LON_PRODUCT", "float")
    add_node_xml(newdoc, root, 'ULY', ULY)
    add_node_xml(newdoc, root, 'ADD_OFFSET', "0")
    add_node_xml(newdoc, root, 'REF_SCALE_FACTOR', "10000")
    add_node_xml(newdoc, root, 'FILLVALUE', "nan")
    add_node_xml(newdoc, root, 'QA_FILLVALUE', "255")
    add_node_xml(newdoc, root, 'MEAN_SUN_AZIMUTH_ANGLE', str(np.mean(solar_azimuth_angle) * 180 / np.pi))
    add_node_xml(newdoc, root, 'MEAN_SUN_ZENITH_ANGLE', str(np.mean(solar_zenithal_angle) * 180 / np.pi))
    add_node_xml(newdoc, root, 'MEAN_VIEW_AZIMUTH_ANGLE', str(np.mean(angle_view_azimuth) * 180 / np.pi))
    add_node_xml(newdoc, root, 'MEAN_VIEW_ZENITH_ANGLE', str(np.mean(angle_view_zenith) * 180 / np.pi))
    add_node_xml(newdoc, root, 'NBAR_SOLAR_ZENITH', str(np.mean(theta_s_out) * 180 / np.pi))
    add_node_xml(newdoc, root, 'ACCODE', "LaSRC version 2.0.1")

    TIRS_SSM_MODEL = find_data_txt(os.path.join(directory_LaSRC, file+"_MTL.txt"), "TIRS_SSM_MODEL", "string")
    add_node_xml(newdoc, root, 'TIRS_SSM_MODEL', TIRS_SSM_MODEL)
    TIRS_SSM_POSITION_STATUS = find_data_txt(os.path.join(directory_LaSRC, file+"_MTL.txt"), "TIRS_SSM_MODEL", "string")
    add_node_xml(newdoc, root, 'TIRS_SSM_POSITION_STATUS', TIRS_SSM_POSITION_STATUS)

    #Save xml file
    f = open(os.path.join(directory_main, file+"xmr.xml"), 'w')
    f.write(newdoc.toxml())
    f.close()


def add_node_xml(newdoc, root, name, value):
    """
    Add node to root element of minidom document newdoc.
    """
    textnode = newdoc.createElement(name)
    root.appendChild(textnode)
    text = newdoc.createTextNode(value)
    textnode.appendChild(text)


def find_data_txt(file, name, type):
    """
    In file, get value in node "name".
    """
    file_txt = open(file)
    for line in file_txt:
        if name in line:
            list_line = str(line).strip().split(" ")
            if type == "string":
                return list_line[2].split('"')[1]
            else:
                return list_line[2]


#Create workspace in input_image
if not os.path.exists(os.path.join(output_image, "workspace")):
    os.makedirs(os.path.join(output_image, "workspace"))

list_image_process = os.listdir(input_initial)
compte = 1
for file in list_image_process:
    print("File {}/{} : {}".format(compte, len(list_image_process), file))
    compte += 1
    t = time.time()
    #Process each image in input_LaSRC_result file
    directory_LaSRC = os.path.join(input_LaSRC_result, file)
    directory_initial = os.path.join(input_initial, file)
    #In workspace, create a file for the image
    if not os.path.exists(os.path.join(output_image, "workspace", file)):
        os.makedirs(os.path.join(output_image, "workspace", file))
    if not os.path.exists(os.path.join(output_image, file)):
        os.makedirs(os.path.join(output_image, file))
    directory_workspace = os.path.join(output_image, "workspace", file)
    directory_main = os.path.join(output_image, file)

    #Open xlsx file
    ws = read_xlsx()
    #Get translation for coregistration
    translation_latitude, translation_longitude = find_translation(file, ws)

    #Calculate angles
    print("Calculating angles...")
    solar_azimuth_angle = open_tif(os.path.join(directory_LaSRC, file + '_solar_azimuth_band4.tif'))
    solar_zenithal_angle = open_tif(os.path.join(directory_LaSRC, file + '_solar_zenith_band4.tif'))
    angle_view_zenith = open_tif(os.path.join(directory_LaSRC, file + '_sensor_zenith_band4.tif'))
    angle_view_azimuth = open_tif(os.path.join(directory_LaSRC, file + '_sensor_azimuth_band4.tif'))
    print("Angles calculated")
    print("")

    #Calculate kernels
    print("Calculating kernels...")
    delta_azimuth = angle_view_azimuth - solar_azimuth_angle
    rows, cols = solar_azimuth_angle.shape
    theta_s_out = find_centre(os.path.join(directory_LaSRC, file + '_MTL.txt'), rows, cols)
    K_geo_sensor = K_geo(angle_view_zenith, solar_zenithal_angle, delta_azimuth)
    K_geo_norm = K_geo(np.zeros((rows, cols)), theta_s_out, np.zeros((rows, cols)))
    K_vol_sensor = K_vol(angle_view_zenith, solar_zenithal_angle, delta_azimuth)
    K_vol_norm = K_vol(np.zeros((rows, cols)), theta_s_out, np.zeros((rows, cols)))
    print("Kernel calculated")
    print("")

    #Process coregistration
    print("Processing coregistration...")
    dst_transform = coregistration(directory_initial, directory_LaSRC, directory_workspace, translation_latitude, translation_longitude, file)
    print("Coregistration processed")
    print("")

    #Save angles
    print("Saving angles...")
    save(os.path.join(directory_main, file + '_solar_azimuth_band4.tif'), solar_azimuth_angle * 180 * 100 / np.pi, dst_transform, np.uint16, 0)
    save(os.path.join(directory_main, file + '_solar_zenith_band4.tif'), solar_zenithal_angle * 180 * 100 / np.pi, dst_transform, np.uint16, 0)
    save(os.path.join(directory_main, file + '_sensor_zenith_band4.tif'), angle_view_zenith * 180 * 100 / np.pi, dst_transform, np.uint16, 0)
    save(os.path.join(directory_main, file + '_sensor_azimuth_band4.tif'), angle_view_azimuth * 180 * 100 / np.pi, dst_transform, np.uint16, 0)
    print("Angles saved")
    print("")

    #Process Fmask
    print("Processing Fmask...")
    process_Fmask(directory_workspace, file, directory_main, dst_transform)
    print("Fmask processed")
    print("")

    #Resampling
    print("Resampling images...")
    resample(directory_workspace, file, directory_main, dst_transform)
    print("Images resampled")
    print("")

    #Process BRDF normalization
    print("Processing BRDF normalization...")
    for band in range(len(list_band_brdf)):
        BRDF(directory_workspace, band, file, K_geo_sensor, K_geo_norm, K_vol_sensor, K_vol_norm, directory_main)
    print("BRDF normalization processed")
    print("")

    print("Processing thermal bands...")
    for band in range(len(list_thermal_band)):
        thermal_band(directory_workspace, directory_main, band, file)
    print("Thermal bands processed")
    print("")

    #Create metadata
    print("Creating metadata...")
    create_xml_metadata(directory_LaSRC, file, directory_main, solar_azimuth_angle, solar_zenithal_angle, angle_view_zenith, angle_view_azimuth, theta_s_out)
    print("Metadata created")
    print("")
    shutil.rmtree(directory_workspace)
    print("Time : ", time.time() - t)
    print("")









