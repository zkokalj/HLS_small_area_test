import os
import rasterio
import numpy as np
from rasterio.windows import Window
from shapely.geometry import box
from rasterio.transform import Affine

#Selection of pixel class used for calculs : with and without water
#comparison_type = "with water"
comparison_type = "without water"

#For images.txt : Sentinel-2 / Landsat 8

#path of folders
input_L8 = "L8"#Landsat 8 images

input_S2 = "S2"#Sentinel-2 images

test = "test"#Results of difference in percentage

#List of bands
list_band = [ ["_sr_band1.tif", "_sr_band1.tif", "A"], ["_sr_band2.tif", "_sr_band2.tif", "A"], ["_sr_band3.tif", "_sr_band3.tif", "A"], ["_sr_band4.tif", "_sr_band4.tif", "A"], ["_sr_band8a.tif", "_sr_band5.tif", "A"], ["_sr_band11.tif", "_sr_band6.tif", "A"], ["_sr_band12.tif", "_sr_band7.tif", "A"], ["_toa_band10.tif", "_toa_band9.tif", "A"]]


def intersection(raster1, raster2):
    """
    Return two objects Rasterio.windows.Window which are intersection of raster1 and raster2 and one object rasterio.transform.Affine
    """
    bb_raster1 = box(raster1.bounds[0], raster1.bounds[1], raster1.bounds[2], raster1.bounds[3])
    bb_raster2 = box(raster2.bounds[0], raster2.bounds[1], raster2.bounds[2], raster2.bounds[3])

    xminR1, yminR1, xmaxR1, ymaxR1 = raster1.bounds
    xminR2, yminR2, xmaxR2, ymaxR2 = raster2.bounds

    intersection = bb_raster1.intersection(bb_raster2)
    transform = Affine(raster1.res[0], 0.0, intersection.bounds[0], 0.0, -raster1.res[1], intersection.bounds[3])

    p1Y = intersection.bounds[3] - raster1.res[1]/2
    p1X = intersection.bounds[0] + raster1.res[0]/2
    p2Y = intersection.bounds[1] + raster1.res[1]/2
    p2X = intersection.bounds[2] - raster1.res[0]/2
    #row index raster1
    row1R1 = int((ymaxR1 - p1Y)/raster1.res[1])
    #row index raster2
    row1R2 = int((ymaxR2 - p1Y)/raster2.res[1])
    #column index raster1
    col1R1 = int((p1X - xminR1)/raster1.res[0])
    #column index raster2
    col1R2 = int((p1X - xminR2)/raster1.res[0])

    #row index raster1
    row2R1 = int((ymaxR1 - p2Y)/raster1.res[1])
    #row index raster2
    row2R2 = int((ymaxR2 - p2Y)/raster2.res[1])
    #column index raster1
    col2R1 = int((p2X - xminR1)/raster1.res[0])
    #column index raster2
    col2R2 = int((p2X - xminR2)/raster1.res[0])

    width1 = col2R1 - col1R1 + 1
    width2 = col2R2 - col1R2 + 1
    height1 = row2R1 - row1R1 + 1
    height2 = row2R2 - row1R2 + 1

    return Window(col1R1, row1R1, width1, height1), Window(col1R2, row1R2, width2, height2), transform

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


def modifier_x(x):
    """
    Used to convert Fmask layer into a filter
    """
    if comparison_type == "without water":
        if x == 0 or x == 16:
            return 1
    else:
        if x == 0 or x == 16 or x == 32 or x == 48:
            return 1
    return 0

#Open images.txt
images_txt = open("images.txt")

number_pixels = [0 for i in range(8)]
difference_value = [0.0 for i in range(8)]
difference_percentage = [0.0 for i in range(8)]

for comparison in images_txt:
    #Get names of image folder
    list_images = comparison.split(";")
    folder_L8 = list_images[1].strip()
    folder_S2 = list_images[0].strip()
    print(folder_S2)
    
    #Open Fmask layers
    Fmask_path_S2 = os.path.join(input_S2, folder_S2, folder_S2+"_Fmask.tif")
    dataset_Fmask_S2 = rasterio.open(Fmask_path_S2)
    Fmask_path_L8 = os.path.join(input_L8, folder_L8, folder_L8+"_Fmask.tif")
    dataset_Fmask_L8 = rasterio.open(Fmask_path_L8)

    #Get intersection between the two images
    emprise_S2, emprise_L8, transform = intersection(dataset_Fmask_S2, dataset_Fmask_L8)
    Fmask_S2 = dataset_Fmask_S2.read(1, window=emprise_S2)
    Fmask_L8 = dataset_Fmask_L8.read(1, window=emprise_L8)

    #Create filter from the two Fmasks 
    vmodifier_x = np.vectorize(modifier_x)
    mask_S2 = vmodifier_x(Fmask_S2)
    mask_L8 = vmodifier_x(Fmask_L8)
    mask = mask_S2 * mask_L8

    for num_band in range(len(list_band)):#For each band
        band = list_band[num_band]
        print(band[0])

        #Open images
        image_path_S2 = os.path.join(input_S2, folder_S2, folder_S2+band[0])
        image_path_L8 = os.path.join(input_L8, folder_L8, folder_L8+band[1])
        dataset_S2 = rasterio.open(image_path_S2)
        dataset_L8 = rasterio.open(image_path_L8)
        image_S2 = dataset_S2.read(1, window=emprise_S2)
        image_L8 = dataset_L8.read(1, window=emprise_L8)

        image_L8 = image_L8 / 10000
        image_S2 = image_S2 / 10000

        
        #Modify spectral band to avoid zero
        masque = np.ones((image_S2.shape))
        masque[image_S2==0]=0
        image_L8_bis = image_L8
        image_S2_bis = image_S2
        image_L8_bis[masque==0]=10
        image_S2_bis[masque==0]=10

        #Calculate difference
        difference = np.abs(image_L8 - image_S2)
        difference_pourcentage = np.abs(image_L8_bis - image_S2_bis) / np.abs(image_S2_bis) * 100

        #Apply masks
        difference_mask =  difference * mask
        difference_pourcentage_masque = difference_pourcentage * mask
        if not os.path.exists(os.path.join(test, folder_S2)):
            os.makedirs(os.path.join(test, folder_S2))
        save(os.path.join(test, folder_S2, band[0]+"_difference_percentage.tif"), difference_pourcentage_masque, transform, np.float32, 0)

        number_pixels[num_band] += np.sum(mask)
        difference_value[num_band] += np.sum(difference_mask)
        difference_percentage[num_band] += np.sum(difference_pourcentage_masque)
        print("Mean difference : ", np.sum(difference_mask) / np.sum(mask))
        print("Mean difference (%) : ", np.sum(difference_pourcentage_masque) / np.sum(mask))
        print("")
        print("")

#Print mean differences 
for i in range(8):
    print(list_band[i][0])
    print("Mean difference on all images : ", difference_value[i] / number_pixels[i])
    print("Mean difference on all images (%) : ", difference_percentage[i] / number_pixels[i])
    print("")