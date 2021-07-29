import os
import rasterio
import numpy as np

#Selection of pixel class used for calculs : with and without water
#comparison_type = "with water"
comparison_type = "without water"

#For images.txt : S30 / Sentinel-2

#path of folders
input_S30 = "S30"#S30 images

input_S2 = "S2"#S2 images

test = "test"#Results of difference in percentage

#List of bands
list_band = [ ["_sr_band1.tif", ".B01.tif", "A"], ["_sr_band2.tif", ".B02.tif", "A"], ["_sr_band3.tif", ".B03.tif", "A"], ["_sr_band4.tif", ".B04.tif", "A"], ["_sr_band5.tif", ".B05.tif", "A"], ["_sr_band6.tif", ".B06.tif", "A"], ["_sr_band7.tif", ".B07.tif", "A"], ["_sr_band8.tif", ".B08.tif", "A"], ["_sr_band8a.tif", ".B8A.tif", "A"], ["_sr_band11.tif", ".B11.tif", "A"], ["_sr_band12.tif", ".B12.tif", "A"], ["_toa_band9.tif", ".B09.tif", "A"], ["_toa_band10.tif", ".B10.tif", "A"], ["_SAA.tif", ".SAA.tif", "B"], ["_SZA.tif", ".SZA.tif", "B"], ["_VAA.tif", ".VAA.tif", "B"], ["_VZA.tif", ".VZA.tif", "B"], ["_Fmask.tif", ".Fmask.tif", "C"]]


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

#Open images.txt
images_txt = open("images.txt")

number_pixels = [0 for i in range(18)]
difference_value = [0.0 for i in range(18)]
difference_percentage = [0.0 for i in range(18)]

for comparison in images_txt:
    #Get names of image folder
    list_images = comparison.split(";")
    folder_S2 = list_images[1].strip()
    folder_S30 = list_images[0].strip()
    
    #Open Fmask layers
    Fmask_path_S30 = os.path.join(input_S30, folder_S30, folder_S30+".Fmask.tif")
    dataset_Fmask_L30 = rasterio.open(Fmask_path_S30)
    Fmask_path_S2 = os.path.join(input_S2, folder_S2, folder_S2+"_Fmask.tif")
    dataset_Fmask_L8 = rasterio.open(Fmask_path_S2)
    Fmask_S30 = dataset_Fmask_L30.read(1)

    #Create filter from the Fmask
    vmodifier_x = np.vectorize(modifier_x)
    mask_S30 = vmodifier_x(Fmask_S30)

    for num_band in range(len(list_band)):#For each band
        band = list_band[num_band]

        #Open images
        image_path_S30 = os.path.join(input_S30, folder_S30, folder_S30+band[1])
        image_path_S2 = os.path.join(input_S2, folder_S2, folder_S2+band[0])
        dataset_S30 = rasterio.open(image_path_S30)
        image_S30 = dataset_S30.read(1)
        dataset_S2 = rasterio.open(image_path_S2)
        image_S2 = dataset_S2.read(1)


        if band[2] == "C":#For Fmask
            image_S30[image_S2==0]=0
            image_S2[image_S30==0]=0
            difference = np.abs(image_S30 - image_S2)
            difference_value[num_band] += np.sum(difference!=0)
            difference_percentage[num_band] += np.sum(difference!=0) * 100
            number_pixels[num_band] += np.sum(mask_S30)
            print("Pour le Fmask, {} % des pixels sont différents.".format(np.sum(difference!=0) / np.sum(mask_S30) * 100))


        if band[2] == "A":#For spectral bands
            image_S30 = image_S30 / 10000
            image_S2 = image_S2 / 10000

            #Modify spectral band to avoid zero
            masque = np.ones((image_S30.shape))
            masque[image_S30==0]=0
            image_S30_bis = image_S30
            image_S2_test_bis = image_S2
            image_S30_bis[masque==0]=10
            image_S2_test_bis[masque==0]=10

            #Calculate difference
            difference = np.abs(image_S30 - image_S2)
            difference_pourcentage = np.abs(image_S30_bis - image_S2_test_bis) / np.abs(image_S30_bis) * 100

            #Apply masks
            difference_mask =  difference * mask_S30
            difference_pourcentage_masque = difference_pourcentage * mask_S30

            #Save difference in percentage in test folder
            if not os.path.exists(os.path.join(test, folder_S2)):
                os.makedirs(os.path.join(test, folder_S2))
            save(os.path.join(test, folder_S2, band[0]+"_difference_percentage.tif"), difference_pourcentage_masque, dataset_S30.transform, np.float32, 0)

            number_pixels[num_band] += np.sum(mask_S30)
            difference_value[num_band] += np.sum(difference_mask)
            difference_percentage[num_band] += np.sum(difference_pourcentage_masque)
            print("Mean difference : ", np.sum(difference_mask) / np.sum(mask_S30))
            print("Mean difference (%) : ", np.sum(difference_pourcentage_masque) / np.sum(mask_S30))
            print("")
            print("")

        if band[2] == "B":#For angle bands
            image_S30[image_S30==40000]=0#In S30, nodata values are 40000
            image_S30[image_S2==0]=0
            image_S2[image_S30==0]=0

            image_S30 = image_S30.astype(np.int32)/100
            image_S2 = image_S2.astype(np.int32)/100

            #Modify spectral band to avoid zero
            masque = np.ones((image_S30.shape))
            masque[image_S30==0]=0
            image_S30_bis = image_S30
            image_S2_bis = image_S2
            image_S30_bis[masque==0]=10
            image_S2_bis[masque==0]=10

            #Calculate difference
            difference = np.abs(image_S30 - image_S2)
            difference_pourcentage = np.abs(image_S30_bis - image_S2_bis) / image_S30_bis * 100

            #Save difference in percentage in test folder
            if not os.path.exists(os.path.join(test, folder_S2)):
                os.makedirs(os.path.join(test, folder_S2))
            save(os.path.join(test, folder_S2, band[0]+"_difference_percentage.tif"), difference_pourcentage, dataset_S30.transform, np.float32, 0)
            
            n, m = difference.shape
            number_pixels[num_band] += n * m
            difference_value[num_band] += np.sum(difference)
            difference_percentage[num_band] += np.sum(difference_pourcentage)
            print("Mean difference : ", np.sum(difference) / (n * m))
            print("Mean difference (%) : ", np.sum(difference_pourcentage) / (n * m))
            print("")
            print("")
            
#Print mean differences  
for i in range(18):
    print(list_band[i][0])
    print("Mean difference on all images : ", difference_value[i] / number_pixels[i])
    print("Mean difference on all images (%) : ", difference_percentage[i] / number_pixels[i])
    print("")