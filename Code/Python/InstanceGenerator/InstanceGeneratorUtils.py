import rasterio
import tifffile as tifffile
from rasterio.mask import mask
import geopandas as gpd
from PIL import Image
import numpy as np
import pandas as pd
import json
import os
import math

import config_instance_generator as config


def clip_tif_by_shp(in_ras, in_shp, out_ras):
    """
    Clips a GeoTIFF raster file using the geometry of a shapefile and saves the resulting cropped image as a new raster.

    Args:
        in_ras (str): Path to the input raster file (GeoTIFF format).
        in_shp (str): Path to the input shapefile used for clipping.
        out_ras (str): Path to save the output cropped raster file.

    Returns:
        None: The resulting cropped raster is saved to the specified output path.
    """
    # Load the shapefile as a GeoDataFrame using GeoPandas
    vector = gpd.read_file(in_shp)

    # Open the input raster file with rasterio
    with rasterio.open(in_ras) as src:
        # Reproject the shapefile to match the coordinate reference system (CRS) of the raster
        vector = vector.to_crs(src.crs)

        # Use the shapefile geometry to mask (clip) the raster
        out_image, out_transform = mask(src, vector.geometry, crop=True)

        # Copy the metadata of the source raster
        out_meta = src.meta.copy()

    # Update the metadata to match the properties of the cropped raster
    out_meta.update({
        "driver": "Gtiff",  # Specify the output format as GeoTIFF
        "height": out_image.shape[1],  # Update the height of the raster
        "width": out_image.shape[2],  # Update the width of the raster
        "transform": out_transform  # Update the affine transform of the raster
    })

    # Write the cropped raster to the output file
    with rasterio.open(out_ras, 'w', **out_meta) as dst:
        dst.write(out_image)


def reshape_matrices_instances(city_folder, city, data_name, data_normalizer, print_log):
    """
    Processes data matrices for a specific city, performing the following steps:
    1. Removes unnecessary dimensions from the raw data files.
    2. Normalizes the data by restoring the correct scale (e.g., temperature values that were scaled by 100).
    3. Removes negative values, setting a minimum value of 0 for all matrices.
    4. Standardizes the size of all matrices to the largest one by scaling smaller matrices.
    5. Exports the processed matrices as TIFF images.

    Args:
        city_folder (str): Path to the folder containing the city's data.
        city (str): Name of the city being processed.
        data_name (list): List of data layer names to process.
        data_normalizer (dict): Dictionary with normalization factors for each data layer.
        print_log (bool): True = log is shown, False = otherwise

    Returns:
        None: The processed matrices are saved as TIFF files in the specified folder.
    """
    # Step 1: Read the raw data layers from TIFF files
    df = {}
    for d in data_name:
        img = rasterio.open(city_folder + '/' + city + '_' + d + config.raw + config.image_format)
        data = img.read()
        df[d] = data

    # Step 2: Remove unnecessary dimensions (keep only the first dimension)
    newDF = {}
    for k in df.keys():
        d = df[k][0]
        newDF[k] = d
    df = newDF

    # Step 3: Normalize the data by multiplying with the normalization factors
    newDF = {}
    for k in df.keys():
        d = np.multiply(df[k], data_normalizer[k], dtype=np.float64)
        newDF[k] = d
    df = newDF

    # Step 4: Remove negative values by setting a minimum of 0
    newDF = {}
    for k in df.keys():
        d = np.maximum(df[k], 0)
        newDF[k] = d
    df = newDF

    # Step 5: Standardize matrix sizes to the largest size
    df = normalize_matrix_size(df)

    # Print information about each processed matrix
    if print_log:
        for k in df.keys():
            print(f"{k} - Shape: {df[k].shape} - MaxVal: {df[k].max()} - MinVal: {df[k].min()}")

    # Step 6: Save the processed matrices as TIFF images
    for k in df.keys():
        dataframe = pd.DataFrame(df[k])  # Convert matrix to DataFrame
        img = Image.fromarray(dataframe.values)  # Convert DataFrame to image
        img.save(city_folder + '/' + city + '_' + k + config.image_format)  # Save as TIFF


def normalize_matrix_size(data):
    """
    Standardizes the size of matrices in the input data by scaling smaller matrices to match the largest one.

    Args:
        data (dict): Dictionary of matrices where the keys are layer names and the values are NumPy arrays.

    Returns:
        dict: Dictionary with all matrices resized to match the largest matrix dimensions.
    """
    # Step 1: Determine the maximum dimensions among all matrices
    maxH, maxW = 0, 0
    for k in data.keys():
        maxH = max(maxH, data[k].shape[0])
        maxW = max(maxW, data[k].shape[1])

    # Step 2: Resize matrices to match the largest dimensions
    new_data = {}
    for k in data.keys():
        if (data[k].shape[0] < maxH) and (data[k].shape[1] < maxW):
            # Create a new matrix with the maximum dimensions
            newDF = np.zeros((maxH, maxW))
            stepH = maxH // data[k].shape[0]  # Scaling factor for height
            stepW = maxW // data[k].shape[1]  # Scaling factor for width

            # Expand the smaller matrix to fit the larger size
            for i in range(data[k].shape[0]):
                for ii in range(stepH):
                    for j in range(data[k].shape[1]):
                        for jj in range(stepW):
                            posH = (i * stepH) + ii
                            posW = (j * stepW) + jj
                            newDF[posH, posW] = data[k][i, j]
            new_data[k] = newDF

        elif data[k].shape[0] > maxH or data[k].shape[1] > maxW:
            print("ERROR: Matrix size exceeds the maximum dimensions.")
        else:
            # If the matrix already matches the largest dimensions, no resizing is needed
            new_data[k] = data[k]

    return new_data


def generate_split_images_by_size(split_size, base_folder, filename, destination_folder, split_label):
    """
    Splits a large image into smaller sub-images (tiles) of approximately equal size based on the specified split size.
    Each tile is saved as a separate image.

    Args:
        split_size (int): Desired approximate size (in pixels) of the smaller sub-images (tiles).
        base_folder (str): Path to the folder containing the original image.
        filename (str): Name of the original image file (including extension).
        destination_folder (str): Path to the folder where the sub-images will be saved.
        split_label (str): Label referring to the split size considered

    Returns:
        None: The sub-images are saved as individual files in the destination folder.
    """
    # Extract the file name (without extension) and the extension of the input image
    name, ext = os.path.splitext(filename)

    # Read the input image using tifffile
    img = tifffile.imread(os.path.join(base_folder, filename))

    # Calculate the number of vertical and horizontal splits based on the split size
    nSplitH = math.ceil(img.shape[0] / split_size)  # Number of splits along the height
    nSplitW = math.ceil(img.shape[1] / split_size)  # Number of splits along the width

    print(f"Size: {split_size} - nSplitH: {nSplitH} - nSplitW: {nSplitW}")

    # Calculate the height and width of each tile
    windowsize_r = int(img.shape[0] / nSplitH) - 1  # Height of each tile
    windowsize_c = int(img.shape[1] / nSplitW) - 1  # Width of each tile

    # Initialize a counter for naming the output tiles
    i = 0

    # Iterate through the image in steps defined by the tile dimensions
    for r in range(0, img.shape[0] - windowsize_r, windowsize_r):  # Loop through rows
        for c in range(0, img.shape[1] - windowsize_c, windowsize_c):  # Loop through columns
            # Extract the sub-image (tile) using slicing
            window = img[r:r + windowsize_r, c:c + windowsize_c]

            # Construct the output file path
            out = os.path.join(destination_folder, f'{name}_{split_label}_{i}{ext}')

            # Save the sub-image (tile) to the destination folder
            tifffile.imwrite(out, window)

            # Increment the tile counter
            i += 1


def create_instances(input_folder, city, data_name, green_type, out_folder, perc_to_be_valid, split_label):
    """
    This method processes a city's folder, extracting all relevant data and information about the impact of
    Nature-Based Solutions (NBS) on different Urban Challenges(UCs).
    For each possible size, it generates an instance and saves it in JSON format.

    Args:
        input_folder (str): Path to the folder containing city data.
        city (str): Name of the city being processed.
        data_name (list): List of data layer names to be read.
        green_type (list): List of NBS types (e.g., Green Roof, Green Wall).
        out_folder (str): Path to the output folder where instances will be saved.
        perc_to_be_valid (float): Minimum percentage of valid positions required for an instance to be considered valid.
        split_label (str): Identifier for the split label (used in file names).

    Returns:
        None: Saves the generated instances directly to the specified output folder.
    """

    # Determine the number of instances based on the presence of urban challenge layers in the input folder
    n_inst = sum(
        config.urban_challenge[0] in file
        for file in os.listdir(input_folder)
        if os.path.isfile(os.path.join(input_folder, file))
    )

    # Iterate through each instance
    df = {}  # Dictionary to hold data layers
    for i in range(0, n_inst, 1):
        print(f"Processing: {city}_{split_label}_{i}")

        # Read each data layer for the current instance
        for d in data_name:
            img = rasterio.open(
                os.path.join(input_folder, f"{city}_{d}_{split_label}_{i}{config.image_format}")
            )
            df[d] = img.read()

        # Check if all required data layers are present and valid
        missing_data = False
        for k in df.keys():
            if np.max(df[k]) <= 0.1:  # Threshold to identify missing data
                missing_data = True
                break

        if not missing_data:
            # Ensure the first dimension is consistent across layers
            newDF = {}
            for k in df.keys():
                d = df[k][0]  # Extract the first dimension of the data
                newDF[k] = d
            df = newDF

            # Create an instance using the data layers and NBS information
            instance = create_instance_from_dict_map_data(df, green_type, perc_to_be_valid)

            # Save the generated instance to the output folder
            if instance is not None:
                instance_path = os.path.join(
                    out_folder, f"{city}_{split_label}_{i}{config.instance_format}"
                )
                with open(instance_path, 'w') as f:
                    json.dump(instance, f)
                print(f"Instance saved: {instance_path}")
        else:
            print("INVALID INSTANCE: Missing data for this instance.")


def create_instance_from_dict_map_data(dataframe, green_type, perc_to_be_valid):
    """
    This method creates an instance based on the provided data and generates a dictionary in JSON format to be saved.

    Args:
        dataframe (dict): A dictionary where keys represent data types and values are corresponding matrices.
        green_type (list): List of Nature-Based Solutions (NBS) types (e.g., Green Roof, Green Wall).
        perc_to_be_valid (float): Minimum percentage of valid positions required for the instance to be considered valid.

    Returns:
        dict or None: Returns the instance dictionary if valid, otherwise None.
    """

    # Extract the first key from the dataframe dictionary, assumed to represent land use data
    k = list(dataframe.keys())[0]  # Key for the land-use matrix
    W = dataframe[k].shape[1]  # Width of the matrix (columns, horizontal, x-axis)
    H = dataframe[k].shape[0]  # Height of the matrix (rows, vertical, y-axis)
    # Note: Matrices in the dataframe are structured as HxW (rows x columns)

    print("W: " + str(W) + " H: " + str(H))

    # Initialize the instance dictionary with metadata and placeholders for forbidden and pre-existing positions
    instance = {
        config.W: W,
        config.H: H,
        config.GT: green_type,
        config.UC: config.urban_challenge,
        config.Forbid: {},
        config.PreExist: {},
        config.A: {}
    }

    # Initialize empty lists for forbidden and pre-existing positions for each NBS type
    for gt in range(len(green_type)):
        instance[config.Forbid][green_type[gt]] = []
        instance[config.PreExist][green_type[gt]] = []

    # Check if the instance has enough valid positions
    land_use_mat = dataframe[k]
    count_forbid = 0

    # Count all forbidden positions
    for i in config.forbidden_for_all:
        result = np.where(land_use_mat == i)
        list_of_coordinates = list(zip(result[0], result[1]))
        count_forbid += len(list_of_coordinates)

    valid_positions = (W * H) - count_forbid  # Calculate the number of valid positions
    if valid_positions > (perc_to_be_valid * W * H):
        # Process land-use data to populate forbidden positions for each green type
        for gt in green_type:
            curr_forb = []

            # Add globally forbidden positions
            for i in config.forbidden_for_all:
                result = np.where(land_use_mat == i)
                list_of_coordinates = list(zip(result[0], result[1]))
                for cord in list_of_coordinates:
                    temp = np.array(cord).tolist()
                    temp = temp[::-1]  # Flip coordinates to x, y format
                    curr_forb.append(temp)

            # Add roof-related forbidden positions (if not GreenRoof)
            if gt != "GreenRoof":
                for i in config.allowed_roof:
                    result = np.where(land_use_mat == i)
                    list_of_coordinates = list(zip(result[0], result[1]))
                    for cord in list_of_coordinates:
                        temp = np.array(cord).tolist()
                        temp = temp[::-1]
                        curr_forb.append(temp)

            # Add wall-related forbidden positions (if not GreenWall)
            if gt != "GreenWall":
                for i in config.allowed_wall:
                    result = np.where(land_use_mat == i)
                    list_of_coordinates = list(zip(result[0], result[1]))
                    for cord in list_of_coordinates:
                        temp = np.array(cord).tolist()
                        temp = temp[::-1]
                        curr_forb.append(temp)

            instance[config.Forbid][gt] = curr_forb

        # Populate pre-existing positions for each green type
        for g in range(len(green_type)):
            curr_pre = []
            for i in config.pre_existent[green_type[g]]:
                result = np.where(land_use_mat == i)
                list_of_coordinates = list(zip(result[0], result[1]))
                for cord in list_of_coordinates:
                    temp = np.array(cord).tolist()
                    temp = temp[::-1]
                    curr_pre.append(temp)
            instance[config.PreExist][green_type[g]] = curr_pre

        # Add urban challenge data to the instance
        for b in config.urban_challenge:
            instance[config.A][b] = dataframe[b].flatten().tolist()  # Flatten data for easier storage in JSON format

        return instance
    else:
        # Print a message if the instance is invalid due to insufficient valid positions
        print("INVALID INSTANCE: " + str(valid_positions) + " valid positions out of " + str(W * H))
        print()
        return None


if __name__ == "__main__":
    # CLIP TIFF WITH SHP COMUNE E CLIP
    for d in config.dataName:
        print("Data: " + d)
        inRas = config.folderBase + config.folderData + d + config.image_format
        for city in config.cities:
            print(city)
            inShp = config.folderBase + config.folderMunicipality + city + config.shp_format
            outRas = config.folderBase + city + '/' + city + '_' + d + config.raw + config.image_format
            clip_tif_by_shp(inRas, inShp, outRas)
        print()

    # NORMALIZE AND RESHAPE THE MATRICES IN ORDER TO HAVE ALL MATRICES OF THE SAME SIZE
    for city in config.cities:
        print(city)
        reshape_matrices_instances(config.folderBase + city, city, config.dataName, config.data_normalizer, True)
        print()

    # SPLIT DATA OF EACH CITY IN TILES
    for city in config.cities:
        print("-----" + city)
        tile_folder = config.folderBase + city + config.folderTiles
        os.makedirs(tile_folder, exist_ok=True)
        for d in config.dataName:
            print("--" + d)
            for i in range(0, len(config.splitSize), 1):
                splitL = config.splitLabel[i]
                outFolder = tile_folder + splitL + "/"
                os.makedirs(outFolder, exist_ok=True)
                generate_split_images_by_size(config.splitSize[i], config.folderBase + city,
                                              city + '_' + d + config.image_format, outFolder, splitL)
        print()

    # GENERATE INSTANCES
    for city in config.cities:
        print(city)
        instance_folder = config.folderBase + city + config.folderInstance
        os.makedirs(instance_folder, exist_ok=True)
        for i in range(0, len(config.splitSize), 1):
            splitL = config.splitLabel[i]
            out_folder = instance_folder + splitL + "/"
            os.makedirs(out_folder, exist_ok=True)
            input_folder = config.folderBase + city + config.folderTiles + splitL + "/"
            create_instances(input_folder, city, config.dataName, config.green_type,
                             out_folder, config.perc_to_be_valid, splitL)
        print()
        print()
