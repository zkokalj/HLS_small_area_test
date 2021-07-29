import os
import rasterio
from rasterio import warp, features
from xml.dom import minidom
import numpy as np
import pyproj
from pyproj.crs import CRS
from pyproj.transformer import Transformer
from scipy import ndimage, interpolate
from scipy.interpolate import lagrange
import numpy.polynomial.polynomial as nppol
from numpy.polynomial.polynomial import Polynomial
import time
import shutil
from rasterio.transform import Affine
import fiona

"""
Process HLS algorithm on Sentinel images. Two directory are mandatory : L1C collection and result of LaSRC algorithm on L1C images.
"""

#Path of L1C images
input_image = "L1C"

#Path of LaSRC result
input_LaSRC_result = "LaSRC_result"

#Path of result images
output_image = "Result"

#Spectral bands processed by BRDF normalization and band adjustment
list_brdf_adjustment = ["_sr_band1.tif", "_sr_band2.tif", "_sr_band3.tif", "_sr_band4.tif", "_sr_band5.tif",  "_sr_band6.tif",  "_sr_band7.tif",  "_sr_band8.tif", "_sr_band8a.tif", "_sr_band11.tif", "_sr_band12.tif"]

#Spectral bands processed only by band adjustment
list_adjustment = ["_sr_band1.tif", "_sr_band2.tif", "_sr_band3.tif", "_sr_band4.tif", "_sr_band8a.tif", "_sr_band11.tif", "_sr_band12.tif"]

#Polynomial coefficients to calculate normalized sun zenithal angle
table_theta_s_out = [31, -0.127, 0.0119, 2.39*10**(-5), -9.46*10**(-7), -1.94*10**(-9), 6.13*10**(-11)]

#Values of BRDF-coefficients
f_iso = [0.0774, 0.0774, 0.1306, 0.1690, 0.2085, 0.2316, 0.2599, 0.3093, 0.3093, 0.3430, 0.2658]
f_geo = [0.0079, 0.0079, 0.0178, 0.0227, 0.0256, 0.0273, 0.0294, 0.0330, 0.0330, 0.0453, 0.0387]
f_vol = [0.0372, 0.0372, 0.0580, 0.0574, 0.0845, 0.1003, 0.1197, 0.1535, 0.1535, 0.1154, 0.0639]

#Band adjustment slope and offset for S2A images
table_value_S2A = [[0.9959, -0.0002], [0.9778, -0.004], [1.0053, -0.0009], [0.9765, 0.0009], [0.9983, -0.0001], [0.9987, -0.0011], [1.003, -0.0012]]

#Band adjustment slope and offset for S2B images
table_value_S2B = [[0.9959, -0.0002], [0.9778, -0.004], [1.0075, -0.0008], [0.9761, 0.001], [0.9966, 0], [1, -0.0003], [0.9867, 0.0004]]

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



def resample_image(directory_LaSRC, directory_workspace, file, directory_input, directory_main):
    """
    Resample each band of an image.
    """
    for name in list_brdf_adjustment:
        dst_transform = resample_band(os.path.join(directory_LaSRC, file+name), os.path.join(directory_workspace, file+name), np.uint16)
        print(os.path.join(directory_LaSRC, file+name))

    #For spectral band 9 and 10
    list_path_initial = os.listdir(os.path.join(directory_input, "GRANULE"))
    path_initial_2 = os.path.join(directory_input, "GRANULE", list_path_initial[0], "IMG_DATA")
    print(os.path.join(directory_LaSRC, file+name))
    for name in os.listdir(path_initial_2):
        if name[-7:]=="B09.jp2":
            dst_transform = resample_band(os.path.join(path_initial_2, name), os.path.join(directory_main, file+"_toa_band9.tif"), np.int16)
        if name[-7:]=="B10.jp2":
            dst_transform = resample_band(os.path.join(path_initial_2, name), os.path.join(directory_main, file+"_toa_band10.tif"), np.int16)
    return dst_transform


def resample_band(input, output, type):
    """
    Resample spectral band.
    """
    dataset = rasterio.open(input)
    dst_transform, width, height = rasterio.warp.aligned_target(dataset.transform, dataset.width, dataset.height, (30, 30))
    destination = np.zeros((height, width), type)
    rasterio.warp.reproject(
            dataset.read(1),
            destination,
            src_transform=dataset.transform,
            src_crs="EPSG:32633",
            dst_transform=dst_transform,
            dst_crs="EPSG:32633",
            resampling=rasterio.warp.Resampling.average)

    #Save spectral band
    with rasterio.open(
                output,
                'w',
                driver='GTiff',
                width=width,
                height=height,
                count=1,
                dtype=type,
                nodata=0,
                transform=dst_transform,
                crs="EPSG:32633") as dst:
        dst.write(destination, indexes=1)
    return dst_transform

def find_solar_angle(path_file, zenith_azimuth):
    """
    Find solar angles (zenithal and azimutal) with information in MTD_TL.xml file. In xml file, grid resolution is 5000 meters.
    """
    #Open MTD_TL.xml file
    xml_data = minidom.parse(os.path.join(path_file, "MTD_TL.xml"))
    root = xml_data.documentElement
    Sun_Angles_Grid = root.getElementsByTagName('Sun_Angles_Grid')
    Zenith = Sun_Angles_Grid[0].getElementsByTagName(zenith_azimuth)
    #Get angle values
    list = []
    for values in Zenith[0].getElementsByTagName('VALUES'):
        list_values = values.childNodes[0].nodeValue.split(" ")
        for i in range(len(list_values)):
            list_values[i] = float(list_values[i])
        list.append(list_values)
    SZA = np.asarray(list)
    #Return a linear interpolation of values. Units is radian.
    return interpolate_angle(SZA) * np.pi / 180
    #return np.ones((3660, 3660)) * np.mean(SZA) * np.pi / 180


def find_view_angle(path_file, directory_input, zenith_azimuth, dst_transform):
    """
    Find view angles (zenithal and azimutal) with information in MTD_TL.xml file. In xml file, grid resolution is 5000 meters.
    For each band, there are five sensors about. Each sensor gives information for only several pixels.
    """

    View_Angles = np.zeros((3660, 3660))

    #Rasterize detector footprint gml
    masque_detecteur = np.int32(rasterize(directory_input, dst_transform))

    #Open MTD_TL.xml
    xml_data = minidom.parse(os.path.join(path_file, "MTD_TL.xml"))
    root = xml_data.documentElement
    Viewing_Incidence_Angles_Grids = root.getElementsByTagName('Viewing_Incidence_Angles_Grids')
    #Keep only angles of spectral band 6
    list_Viewing_Incidence = [i for i in Viewing_Incidence_Angles_Grids if i.getAttribute("bandId")=="5"]

    #For each detector, view angles are calculated
    for detector in list_Viewing_Incidence:
        #Get values for detector
        detector_id = int(detector.getAttribute("detectorId"))
        Zenith = detector.getElementsByTagName(zenith_azimuth)
        list = []
        for values in Zenith[0].getElementsByTagName('VALUES'):
            list_values = values.childNodes[0].nodeValue.split(" ")
            for i in range(len(list_values)):
                list_values[i] = float(list_values[i])
            list.append(list_values)
        detector_angle = np.asarray(list)
        x = np.arange(0, 23, 1)
        y = np.arange(0, 23, 1)
        xx, yy = np.meshgrid(x, y)
        #Add columns to avoid nan values in interpolation
        vadd_column = np.vectorize(add_column)
        vadd_column.excluded.add(2)
        detector_angle = vadd_column(yy, xx, detector_angle)

        #Interpolate values
        f = interpolate.interp2d(x, y, detector_angle, kind='linear')
        x_new = np.linspace(0, 22, 3660)
        y_new = np.linspace(0, 22, 3660)
        array_interpolate = f(x_new, y_new)

        #Calculate mask to separate value dealt with by this detector from other values.
        specific_mask_detector = np.copy(masque_detecteur)

        specific_mask_detector.astype(np.int32)
        for i in range(1, 13):
            if i != detector_id:
                specific_mask_detector = (masque_detecteur - i) * specific_mask_detector
        if np.max(np.abs(specific_mask_detector)) > 0:
            specific_mask_detector = np.abs(specific_mask_detector / (np.max(np.abs(specific_mask_detector))))
            #Apply mask on interpolated values
            detector_angle_interpolated = specific_mask_detector * array_interpolate
            View_Angles = View_Angles + detector_angle_interpolated
    #Return a linear interpolation of values. Units is radian
    return View_Angles * np.pi / 180


def rasterize(directory_input, dst_transform):
    """
    Rasterize detector footprint.
    """
    temp = os.listdir(os.path.join(directory_input, "GRANULE"))
    with fiona.open(os.path.join(directory_input, "GRANULE", temp[0], "QI_DATA", "MSK_DETFOO_B06.gml"), "r") as shapefile:
        shapes = [(feature["geometry"], int(feature["properties"]["gml_id"][-4:-2])) for feature in shapefile]
        image = features.rasterize(
            ((g, v) for g, v in shapes),
            out_shape=(3660, 3660),
            transform=dst_transform)
    return image


def add_column(xx, yy, detector_angle):
    """
    To avoid nan values in interpolation, values are added in each side of detector.
    """
    if np.isnan(detector_angle[xx][yy]):
        if yy < 22:
            if not np.isnan(detector_angle[xx][yy+1]):
                return detector_angle[xx][yy+1]
        if yy > 0:
            if not np.isnan(detector_angle[xx][yy-1]):
                return detector_angle[xx][yy-1]
        return 0.0
    else:
        return detector_angle[xx][yy]

def interpolate_angle(array):
    """
    Interpolate a 3660 * 3660 grid from a 23 * 23 grid with linear interpolation
    """
    x = np.arange(0, 23, 1)
    y = np.arange(0, 23, 1)
    f = interpolate.interp2d(x, y, array, kind='linear', fill_value=np.nan)
    x_new = np.linspace(0, 22, 3660)
    y_new = np.linspace(0, 22, 3660)
    array_interpolate = f(x_new, y_new)
    return array_interpolate

def save_angles(solar_azimuth_angle, solar_zenithal_angle, angle_view_zenith, angle_view_azimuth, dst_transform, directory_main, file):
    """
    Save angles values. Values are saved in degrees.
    """
    save(os.path.join(directory_main, file+"_SAA.tif"), solar_azimuth_angle * 100 * 180 / np.pi, dst_transform, np.uint16, 0)
    save(os.path.join(directory_main, file+"_SZA.tif"), solar_zenithal_angle * 100 * 180 / np.pi, dst_transform, np.uint16, 0)
    save(os.path.join(directory_main, file+"_VZA.tif"), angle_view_zenith * 100 * 180 / np.pi, dst_transform, np.uint16, 0)
    save(os.path.join(directory_main, file+"_VAA.tif"), angle_view_azimuth * 100 * 180 / np.pi, dst_transform, np.uint16, 0)

def find_centre(path_file):
    """
    Find latitude of center of the image with information in MTD_TL.xml file.
    """
    #Open MTD_TL.xml file
    xml_data = minidom.parse(os.path.join(path_file, "MTD_TL.xml"))
    root = xml_data.documentElement
    Sun_Angles_Grid_Y = root.getElementsByTagName('ULY')
    top_Y = Sun_Angles_Grid_Y[0].childNodes[0].nodeValue
    center_Y = int(top_Y) - 30 * 3660 / 2

    Sun_Angles_Grid_X = root.getElementsByTagName('ULX')
    top_X = Sun_Angles_Grid_X[0].childNodes[0].nodeValue
    center_X = int(top_X) + 30 * 3660 / 2

    #Transform projection into latitude/longitude
    crs = CRS(proj="longlat", ellps="WGS84")
    transformer = Transformer.from_crs("epsg:32633", crs)
    longitude, latitude = transformer.transform(center_X, center_Y)

    #Calculate polynom
    sum = 0
    for i in range(len(table_theta_s_out)):
        sum += table_theta_s_out[i] * latitude**i
    return np.ones((3660, 3660)) * sum * np.pi / 180

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
        print("Alert cos_t !!!")


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

def adjacent_cloud_Sentinel(Fmask, cloud_shadow, dilatation):
    """
    Return Fmask with all classes (adjacent cloud included).
    """
    cloud_shadow_Fmask = cloud_shadow * Fmask
    dilatation_Fmask = (dilatation - cloud_shadow) * 4
    other_Mask = - (dilatation - 1)
    remain_Fmask = other_Mask * Fmask
    return cloud_shadow_Fmask + dilatation_Fmask + remain_Fmask

def convert_bit_index(x):
    """
    Convert x in a bit index.
    """
    if x == 666666666:#if x is a non data value
        return 255
    x_string = str(x)
    sum = 0
    for i in range(1,6):
        if str(i) in x_string:
            sum += 2**i
    return sum

def process_Fmask(directory_LaSRC, directory_main, file, dst_transform):
    """
    Process Fmask : modify values, add adjacent cloud classes, transform values into a bit index.
    """
    path_Fmask = os.path.join(directory_LaSRC, file+"_Fmask4.tif")
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
        dilatation = ndimage.binary_dilation(cloud_shadow, structure=structure, iterations=15).astype(cloud_shadow.dtype)

        #Modify value after dilatation adding class adjacent cloud / cloud shadow
        Fmask2 = adjacent_cloud_Sentinel(Fmask, cloud_shadow, dilatation)

        #Resample Fmask with a cubic convolution adapted for Fmask
        channel_small_values = nppol.polyval(Fmask2, small_values)
        channel_small_values = (channel_small_values+0.1).astype("i8")
        filtre = np.array([[10**8, 10**7, 10**6],[10**5, 10**4, 10**3],[10**2, 10, 1]])
        resample = ndimage.convolve(channel_small_values, filtre)
        resample1 = resample[1:, 1:]
        resampled_image = resample1[::3, ::3]
        vFmask_create_bit_index = np.vectorize(convert_bit_index)
        Fmask_bit_index = vFmask_create_bit_index(resampled_image)

        #save image
        save(os.path.join(directory_main, file+"_Fmask.tif"), Fmask_bit_index, dst_transform, np.int16, 255)

def BRDF(path_file, num_band, file, K_geo_sensor, K_geo_norm, K_vol_sensor, K_vol_norm):
    """
    Process BRDF-normalization on band num_band.
    """
    path_band = os.path.join(path_file, file+list_brdf_adjustment[num_band])
    if os.path.exists(path_band):
        channel = rasterio.open(path_band).read(1)
        nan_array = channel / 10000
        nan_array[channel == 0] = np.nan

        f_iso_l = f_iso[num_band]
        f_geo_l = f_geo[num_band]
        f_vol_l = f_vol[num_band]

        #Calculate c-factor
        c = ( f_iso_l + f_geo_l * K_geo_norm + f_vol_l * K_vol_norm) / ( f_iso_l + f_geo_l * K_geo_sensor + f_vol_l * K_vol_sensor )

        #Calculate BRDF-normalized image
        BRDF_image = c * nan_array
        return BRDF_image

def band_adjustment(file, array, num_band):
    """
    Calculate band adjustment for each band in file.
    """
    #Values depends on S2A or S2B images
    if file[2] == "A":
        table = table_value_S2A
    else:
        table = table_value_S2B

    adjustment = array * table[num_band][0] + table[num_band][1]
    return adjustment


def create_xml_metadata(directory_LaSRC, solar_azimuth_angle, solar_zenithal_angle, angle_view_zenith, angle_view_azimuth, theta_s_out, file, directory_main):
    """
    Create xml file with metadata
    """
    #Open LaSRC metadata
    xml_data_MTD_MSIL1C = minidom.parse(os.path.join(directory_LaSRC, "MTD_MSIL1C.xml"))
    xml_data_MTD_TL = minidom.parse(os.path.join(directory_LaSRC, "MTD_TL.xml"))

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
    PRODUCT_URI = xml_data_MTD_MSIL1C.getElementsByTagName('PRODUCT_URI')
    add_node_xml(newdoc, root, 'PRODUCT_URI', PRODUCT_URI[0].firstChild.nodeValue)
    SENSING_TIME = xml_data_MTD_MSIL1C.getElementsByTagName('PRODUCT_START_TIME')
    add_node_xml(newdoc, root, 'SENSING_TIME', SENSING_TIME[0].firstChild.nodeValue)
    add_node_xml(newdoc, root, 'SPATIAL_COVERAGE', str(100 - count_Fmask_None / (n*m) * 100))
    add_node_xml(newdoc, root, 'CLOUD_COVERAGE', str((count_Fmask_cloud + count_Fmask_cloud_shadow) / (n*m - count_Fmask_None) * 100))
    ULX = xml_data_MTD_TL.getElementsByTagName('ULX')
    add_node_xml(newdoc, root, 'ULX', ULX[0].firstChild.nodeValue)
    ULY = xml_data_MTD_TL.getElementsByTagName('ULY')
    add_node_xml(newdoc, root, 'ULY', ULY[0].firstChild.nodeValue)
    add_node_xml(newdoc, root, 'SPATIAL_RESAMPLING_ALG', "area weighted average")
    add_node_xml(newdoc, root, 'ADD_OFFSET', "0")
    add_node_xml(newdoc, root, 'REF_SCALE_FACTOR', "10000")
    add_node_xml(newdoc, root, 'FILLVALUE', "0")
    add_node_xml(newdoc, root, 'QA_FILLVALUE', "255")
    add_node_xml(newdoc, root, 'MEAN_SUN_AZIMUTH_ANGLE', str(np.mean(solar_azimuth_angle) * 180 / np.pi))
    add_node_xml(newdoc, root, 'MEAN_SUN_ZENITH_ANGLE', str(np.mean(solar_zenithal_angle) * 180 / np.pi))
    add_node_xml(newdoc, root, 'MEAN_VIEW_AZIMUTH_ANGLE', str(np.mean(angle_view_azimuth) * 180 / np.pi))
    add_node_xml(newdoc, root, 'MEAN_VIEW_ZENITH_ANGLE', str(np.mean(angle_view_zenith) * 180 / np.pi))
    add_node_xml(newdoc, root, 'NBAR_SOLAR_ZENITH', str(np.mean(theta_s_out) * 180 / np.pi))

    list_spectral_band = ["01", "02", "03", "04", "8A", "11", "12"]
    if file[2] == "A":
        table = table_value_S2A
    else:
        table = table_value_S2B
    for i in range(7):
        add_node_xml(newdoc, root, "MSI_BAND_{}_BANDPASS_ADJUSTMENT_SLOPE".format(list_spectral_band[i]), str(table[i][0]))
        add_node_xml(newdoc, root, "MSI_BAND_{}_BANDPASS_ADJUSTMENT_OFFSET".format(list_spectral_band[i]), str(table[i][1]))
    add_node_xml(newdoc, root, 'ACCODE', "LaSRC version 2.0.1")

    #Save xml file
    f = open(os.path.join(directory_main, file+"_xmr.xml"), 'w')
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

#Create workspace in input_image
if not os.path.exists(os.path.join(output_image, "workspace")):
    os.makedirs(os.path.join(output_image, "workspace"))

list_image_process = os.listdir(input_LaSRC_result)
compte = 1
for file in list_image_process:
    print("File {}/{} : {}".format(compte, len(list_image_process), file))
    compte += 1
    t = time.time()
    #Process each image in input_LaSRC_result file
    directory_input = os.path.join(input_image, file+".SAFE")
    directory_LaSRC = os.path.join(input_LaSRC_result, file)
    #In workspace, create a file for the image
    if not os.path.exists(os.path.join(output_image, "workspace", file)):
        os.makedirs(os.path.join(output_image, "workspace", file))
    if not os.path.exists(os.path.join(output_image, file)):
        os.makedirs(os.path.join(output_image, file))
    directory_workspace = os.path.join(output_image, "workspace", file)
    directory_main = os.path.join(output_image, file)

    #Resample image
    print("Resampling images...")
    dst_transform = resample_image(directory_LaSRC, directory_workspace, file, directory_input, directory_main)
    print("Images resampled")
    print("")

    #Calculate values common to each band of an image
    print("Angles calculating...")
    solar_azimuth_angle = find_solar_angle(directory_LaSRC, 'Azimuth')
    solar_zenithal_angle = find_solar_angle(directory_LaSRC, 'Zenith')
    angle_view_zenith = find_view_angle(directory_LaSRC, directory_input, 'Zenith', dst_transform)
    angle_view_azimuth = find_view_angle(directory_LaSRC, directory_input, 'Azimuth', dst_transform)
    print("Angles calculated")
    print("")

    #Save angles
    print("Saving angles...")
    save_angles(solar_azimuth_angle, solar_zenithal_angle, angle_view_zenith, angle_view_azimuth, dst_transform, directory_main, file)
    print("Angles saved")
    print("")

    #Calculate kernels
    print("Calculating kernels...")
    delta_azimuth = angle_view_azimuth - solar_azimuth_angle
    theta_s_out = find_centre(directory_LaSRC)
    K_geo_sensor = K_geo(angle_view_zenith, solar_zenithal_angle, delta_azimuth)
    K_geo_norm = K_geo(np.zeros((3660, 3660)), theta_s_out, np.zeros((3660, 3660)))
    K_vol_sensor = K_vol(angle_view_zenith, solar_zenithal_angle, delta_azimuth)
    K_vol_norm = K_vol(np.zeros((3660, 3660)), theta_s_out, np.zeros((3660, 3660)))
    print("Kernels calculated")
    print("")

    #Process Fmask
    print("Processing Fmask...")
    process_Fmask(directory_LaSRC, directory_main, file, dst_transform)
    print("Fmask processed")
    print("")

    #Create metadata
    print("Creating metadata...")
    create_xml_metadata(directory_LaSRC, solar_azimuth_angle, solar_zenithal_angle, angle_view_zenith, angle_view_azimuth, theta_s_out, file, directory_main)
    print("Metadata created...")
    print("")

    #Process BRDF normalization and band adjustment
    print("Processing BRDF normalization and band adjustment...")
    for num_band in range(len(list_brdf_adjustment)):
        array = BRDF(directory_workspace, num_band, file, K_geo_sensor, K_geo_norm, K_vol_sensor, K_vol_norm)
        if list_brdf_adjustment[num_band] in list_adjustment:
            index = list_adjustment.index(list_brdf_adjustment[num_band])
            array = band_adjustment(file, array, index)
        array = array * 10000
        array[array==np.nan]=-9999
        save(os.path.join(directory_main, file+list_brdf_adjustment[num_band]), array, dst_transform, np.int16, -9999)
    print("BRDF normalization and band adjustment processed...")

    shutil.rmtree(directory_workspace)
    print("Time : ", time.time() - t)
    print("")





