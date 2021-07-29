import os
import rasterio
import numpy as np
from rasterio.windows import Window
from shapely.geometry import box
from rasterio.transform import Affine

#Selection of pixel class used for calculs : with and without water
#comparison_type = "with water"
comparison_type = "without water"


#For images.txt : Landsat 5 / Landsat 7

#path of folders
input_L5 = "L5"#Landsat 5 images

input_L7 = "L7"#Landsat 7 images

test = "test"#Results of difference in percentage

#List of bands
list_band = ["_sr_band1.tif", "_sr_band2.tif", "_sr_band3.tif", "_sr_band4.tif", "_sr_band5.tif", "_bt_band6.tif", "_sr_band7.tif"]



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

    p1Y = intersection.bounds[3] - raster1.res[1]/2 -30000
    p1X = intersection.bounds[0] + raster1.res[0]/2 + 30000
    p2Y = intersection.bounds[1] + raster1.res[1]/2 + 80000
    p2X = intersection.bounds[2] - raster1.res[0]/2 - 30000
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

number_pixels = [0 for i in range(7)]
difference_value = [0.0 for i in range(7)]
difference_percentage = [0.0 for i in range(7)]

for comparison in images_txt:
    #Get names of image folder
    list_images = comparison.split(";")
    folder_L7 = list_images[1].strip()
    folder_L5 = list_images[0].strip()
    print(folder_L5)
    
    #Open Fmask layers
    Fmask_path_L5 = os.path.join(input_L5, folder_L5, folder_L5+"_Fmask.tif")
    dataset_Fmask_L5 = rasterio.open(Fmask_path_L5)
    Fmask_path_L7 = os.path.join(input_L7, folder_L7, folder_L7+"_Fmask.tif")
    dataset_Fmask_L7 = rasterio.open(Fmask_path_L7)

    #Get intersection between the two images
    emprise_L5, emprise_L7, transform = intersection(dataset_Fmask_L5, dataset_Fmask_L7)
    Fmask_L5 = dataset_Fmask_L5.read(1, window=emprise_L5)
    Fmask_L7 = dataset_Fmask_L7.read(1, window=emprise_L7)

    #Create filter from the two Fmasks 
    vmodifier_x = np.vectorize(modifier_x)
    mask_L5 = vmodifier_x(Fmask_L5)
    mask_L7 = vmodifier_x(Fmask_L7)
    mask = mask_L5 * mask_L7

    for num_band in range(len(list_band)):#For each band
        print(list_band[num_band])

        #Open images
        image_path_L5 = os.path.join(input_L5, folder_L5, folder_L5+list_band[num_band])
        image_path_L7 = os.path.join(input_L7, folder_L7, folder_L7+list_band[num_band])
        dataset_L5 = rasterio.open(image_path_L5)
        dataset_L7 = rasterio.open(image_path_L7)
        image_L5 = dataset_L5.read(1, window=emprise_L5)
        image_L7 = dataset_L7.read(1, window=emprise_L7)

        #Modify spectral band to avoid zero
        masque = np.ones((image_L5.shape))
        masque[image_L5==0]=0
        masque[image_L7==0]=0
        image_L5_bis = image_L5
        image_L7_test_bis = image_L7
        image_L5_bis[masque==0]=10
        image_L7_test_bis[masque==0]=10

        #Calculate difference
        difference = np.abs(image_L5 - image_L7)
        difference_pourcentage = np.abs(image_L5_bis - image_L7_test_bis) / np.abs(image_L5_bis) * 100

        #Apply masks
        difference_mask =  difference * masque * mask
        difference_pourcentage_masque = difference_pourcentage * masque * mask

        #Save difference in percentage in test folder
        if not os.path.exists(os.path.join(test, folder_L5)):
                os.makedirs(os.path.join(test, folder_L5))
        save(os.path.join(test, folder_L5, list_band[num_band]+"_difference_percentage_masque.tif"), difference_pourcentage_masque, transform, np.float32, 0)

        number_pixels[num_band] += np.sum(mask)
        difference_value[num_band] += np.sum(difference_mask)
        difference_percentage[num_band] += np.sum(difference_pourcentage_masque)
        print("Mean difference : ", np.sum(difference_mask) / np.sum(masque))
        print("Mean difference (%) : ", np.sum(difference_pourcentage_masque) / np.sum(masque))
        print("")
        print("")
        
#Print mean differences       
for i in range(7):
    print(list_band[i])
    print("Mean difference on all images : ", difference_value[i] / number_pixels[i])
    print("Mean difference on all images (%) : ", difference_percentage[i] / number_pixels[i])
    print("")