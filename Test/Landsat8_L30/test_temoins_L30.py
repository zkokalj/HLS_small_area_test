import os
import rasterio
import numpy as np
from rasterio.windows import Window
from shapely.geometry import box
from rasterio.transform import Affine

#Selection of pixel class used for calculs : with and without water
#comparison_type = "with water"
comparison_type = "without water"

#For images.txt : L30 / Landsat 8

#path of folders
input_L30 = "L30"#L30 images

input_L8 = "L8"#Landsat 8 images

test = "test"#Results of difference in percentage

list_band = [ ["_sr_band1.tif", ".B01.tif", "A"], ["_sr_band2.tif", ".B02.tif", "A"], ["_sr_band3.tif", ".B03.tif", "A"], ["_sr_band4.tif", ".B04.tif", "A"], ["_sr_band5.tif", ".B05.tif", "A"], ["_sr_band6.tif", ".B06.tif", "A"], ["_sr_band7.tif", ".B07.tif", "A"], ["_bt_band11.tif", ".B11.tif", "A"], ["_toa_band9.tif", ".B09.tif", "A"], ["_bt_band10.tif", ".B10.tif", "A"], ["_solar_azimuth_band4.tif", ".SAA.tif", "B"], ["_solar_zenith_band4.tif", ".SZA.tif", "B"], ["_sensor_azimuth_band4.tif", ".VAA.tif", "B"], ["_sensor_zenith_band4.tif", ".VZA.tif", "B"], ["_Fmask.tif", ".Fmask.tif", "C"]]

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
        if (x == 0 or x == 16 ):
            return 1
    else:
        if x == 0 or x == 16 or x == 32 or x == 48:
            return 1
    return 0

#Open images.txt
images_txt = open("images.txt")

number_pixels = [0 for i in range(15)]
difference_value = [0.0 for i in range(15)]
difference_percentage = [0.0 for i in range(15)]

for comparison in images_txt:
    #Get names of image folder
    list_images = comparison.split(";")
    folder_L8 = list_images[1].strip()
    folder_L30 = list_images[0].strip()
    
    #Open Fmask layers
    Fmask_path_L30 = os.path.join(input_L30, folder_L30, folder_L30+".Fmask.tif")
    dataset_Fmask_L30 = rasterio.open(Fmask_path_L30)
    Fmask_path_L8 = os.path.join(input_L8, folder_L8, folder_L8+"_Fmask.tif")
    dataset_Fmask_L8 = rasterio.open(Fmask_path_L8)

    #Get intersection between the two images
    emprise_L30, emprise_L8, transform = intersection(dataset_Fmask_L30, dataset_Fmask_L8)
    Fmask_L30 = dataset_Fmask_L30.read(1, window=emprise_L30)
    
    #Create filter from the Fmask
    vmodifier_x = np.vectorize(modifier_x)
    mask_L30 = vmodifier_x(Fmask_L30)

    for num_band in range(len(list_band)):#For each band
        band = list_band[num_band]
        print(list_band[num_band][0])

        #Open images
        image_path_L30 = os.path.join(input_L30, folder_L30, folder_L30+list_band[num_band][1])
        image_path_L8 = os.path.join(input_L8, folder_L8, folder_L8+list_band[num_band][0])
        dataset_L30 = rasterio.open(image_path_L30)
        dataset_L8 = rasterio.open(image_path_L8)
        image_L30 = dataset_L30.read(1, window=emprise_L30)
        image_L8 = dataset_L8.read(1, window=emprise_L8)

        if band[2] == "C":#For Fmask
            image_L30[image_L8==0]=0
            image_L8[image_L30==0]=0
            difference = np.abs(image_L30 - image_L8)
            difference_value[num_band] += np.sum(difference!=0)
            difference_percentage[num_band] += np.sum(difference!=0) * 100
            number_pixels[num_band] += np.sum(mask_L30)
            print("Pour le Fmask, {} % des pixels sont différents.".format(np.sum(difference!=0) / np.sum(mask_L30) * 100))
        
        if band[2] == "A":#For spectral bands
            image_L30[image_L30==-9999]=0#In L30, nodata value is -9999
            image_L30[image_L8==0]=0
            image_L8[image_L30==0]=0
            if band[1] == "B10.tif" or band[1] == "B11.tif":
                image_L30 = image_L30 / 100
                image_L8 = image_L8 / 100
            else:
                image_L30 = image_L30 / 10000
                image_L8 = image_L8 / 10000

            
            #Modify spectral band to avoid zero
            masque = np.ones((image_L30.shape))
            masque[image_L30==0]=0
            image_L30_bis = image_L30
            image_L8_test_bis = image_L8
            image_L30_bis[masque==0]=10
            image_L8_test_bis[masque==0]=10

            #Calculate difference
            difference = np.abs(image_L30 - image_L8)
            difference_pourcentage = np.abs(image_L30_bis - image_L8_test_bis) / np.abs(image_L30_bis) * 100

            #Apply masks
            difference_mask =  difference * mask_L30
            difference_pourcentage_masque = difference_pourcentage * mask_L30

            #Save difference in percentage in test folder
            if not os.path.exists(os.path.join(test, folder_L8)):
                os.makedirs(os.path.join(test, folder_L8))
            save(os.path.join(test, folder_L8, band[0]+"_difference_percentage.tif"), difference_pourcentage_masque, transform, np.float32, 0)

            number_pixels[num_band] += np.sum(mask_L30)
            difference_value[num_band] += np.sum(difference_mask)
            difference_percentage[num_band] += np.sum(difference_pourcentage_masque)
            print("Mean difference : ", np.sum(difference_mask) / np.sum(mask_L30))
            print("Mean difference (%) : ", np.sum(difference_pourcentage_masque) / np.sum(mask_L30))
            print("")
            print("")

        if band[2] == "B":#For angle bands
            image_L30[image_L30==40000]=0#In S30, nodata values are 40000
            image_L30[image_L8==0]=0
            image_L8[image_L30==0]=0

            image_L30 = image_L30.astype(np.int32)/100
            image_L8 = image_L8.astype(np.int32)/100

            #Modify spectral band to avoid zero
            masque = np.ones((image_L30.shape))
            masque[image_L30==0]=0
            image_L30_bis = image_L30
            image_L8_bis = image_L8
            image_L30_bis[masque==0]=10
            image_L8_bis[masque==0]=10

            #Calculate difference
            difference = np.abs(image_L30 - image_L8)
            difference_pourcentage = np.abs(image_L30_bis - image_L8_bis) / image_L30_bis * 100

            #Save difference in percentage in test folder
            if not os.path.exists(os.path.join(test, folder_L8)):
                os.makedirs(os.path.join(test, folder_L8))
            save(os.path.join(test, folder_L8, band[0]+"_difference_percentage.tif"), difference_pourcentage, transform, np.float32, 0)
            
            n, m = difference.shape
            number_pixels[num_band] += n * m
            difference_value[num_band] += np.sum(difference)
            difference_percentage[num_band] += np.sum(difference_pourcentage)
            print("Mean difference : ", np.sum(difference) / (n * m))
            print("Mean difference (%) : ", np.sum(difference_pourcentage) / (n * m))
            print("")
            print("")
            
#Print mean differences  
for i in range(15):
    print(list_band[i][0])
    print("Mean difference on all images : ", difference_value[i] / number_pixels[i])
    print("Mean difference on all images (%) : ", difference_percentage[i] / number_pixels[i])
    print("")

































