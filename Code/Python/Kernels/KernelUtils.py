import numpy as np
import json

max_val = 1.0  # kernel default max value
min_val = 0.1  # kernel default min value


def generate_kernel_2d(w, h, generation_alg, min_val, max_val):
    """
    Redirect the generation of a 2D kernel based on the specified generation algorithm:
        0: generate_2d_all_one,  # Kernel with all values set to 1
        1: generate_2d_all_zero,  # Kernel with all values set to 0
        2: generate_2d_random,  # Kernel with random values
        3: generate_2d_linear_from_max_to_min,  # Linearly decreasing kernel from a central maximum value
        4: generate_2d_linear_from_min_to_max,  # Linearly increasing kernel from a central minimum value
        5: generate_2d_point_effect,  # Kernel with a central peak

    Args:
        w (int): Width of the kernel.
        h (int): Height of the kernel.
        generation_alg (int): Index specifying the generation method.
        min_val (float): Minimum value used in some algorithms.
        max_val (float): Maximum value used in some algorithms.

    Returns:
        np.ndarray: The generated 2D kernel.
    """
    # Dictionary mapping algorithm indices to their corresponding functions.
    switcher = {
        0: generate_2d_all_one,  # Kernel with all values set to 1
        1: generate_2d_all_zero,  # Kernel with all values set to 0
        2: generate_2d_random,  # Kernel with random values
        3: generate_2d_linear_from_max_to_min,  # Linearly decreasing kernel
        4: generate_2d_linear_from_min_to_max,  # Linearly increasing kernel
        5: generate_2d_point_effect,  # Kernel with a central peak
    }
    # Retrieve the function corresponding to the chosen algorithm,
    # defaulting to the all-one kernel.
    func = switcher.get(generation_alg, generate_2d_all_one)
    # Call the selected function with the provided parameters.
    return func(w, h, min_val, max_val)


def generate_2d_all_one(w, h):
    """
    Generates a 2D kernel of size w x h filled with the value 1.
    """
    return np.ones((w, h), dtype=float)


def generate_2d_all_zero(w, h):
    """
    Generates a 2D kernel of size w x h filled with the value 0.
    """
    return np.zeros((w, h), dtype=float)


def generate_2d_random(w, h, min_val, max_val):
    """
    Generates a 2D kernel of size w x h with random values
    between min_val and max_val.
    """
    return np.random.uniform(min_val, max_val, (w, h))


def generate_2d_point_effect(w, h):
    """
    Generates a 2D kernel with a central peak value of 1 and
    all other values set to 0.
    """
    k = np.zeros((w, h), dtype=float)
    mid_row = w // 2  # Calculate the middle row
    mid_col = h // 2  # Calculate the middle column
    k[mid_row][mid_col] = 1  # Set the central peak
    return k


def generate_2d_linear_from_max_to_min(w, h, min_val, max_val):
    """
    Generates a 2D kernel where values decrease linearly
    from the center (max_val) to the edges (min_val).
    """
    k = np.zeros((w, h), dtype=float)
    mid_row, mid_col = w // 2, h // 2  # Center of the kernel
    max_dist = max(mid_row, mid_col)  # Maximum distance from the center
    for i in range(w):
        for j in range(h):
            dist = max(abs(i - mid_row), abs(j - mid_col))  # Current distance
            k[i][j] = max_val - (max_val - min_val) * dist / max_dist  # Compute value
    return k


def generate_2d_linear_from_min_to_max(w, h, min_val, max_val):
    """
    Generates a 2D kernel where values increase linearly
    from the center (min_val) to the edges (max_val).
    """
    k = np.zeros((w, h), dtype=float)
    mid_row, mid_col = w // 2, h // 2  # Center of the kernel
    max_dist = max(mid_row, mid_col)  # Maximum distance from the center
    for i in range(w):
        for j in range(h):
            dist = max(abs(i - mid_row), abs(j - mid_col))  # Current distance
            k[i][j] = min_val + (max_val - min_val) * dist / max_dist  # Compute value
    return k


def generate_kernels_2d(alg, green_types, urban_challenges, sizes, ranges):
    """
    Generates 2D kernels for each green types, urban challenges given the relative sizes and ranges of values.

    Args:
        alg (int): Algorithm index for kernel generation.
        green_types (list): List of green type identifiers.
        urban_challenges (list): List of urban challenge type identifiers.
        sizes (dict): Dictionary where keys are green_types and values are (rows, cols) tuples.
        ranges (dict): Dictionary where keys are green_types and values are (min_val, max_val) tuples.

    Returns:
        dict: A dictionary containing the generated kernels and their metadata.
    """
    kernel = {}

    # Save the list of green types
    kernel["GreenType"] = list(green_types)

    # Save the list of urban challenge
    kernel["UrbanChallenge"] = list(urban_challenges)

    # Save the sizes of the kernels
    kernel["SizeK"] = {}
    for b_type in urban_challenges:
        kernel["SizeK"][b_type] = {}
        for g_type in green_types:
            kernel["SizeK"][b_type][g_type] = {
                "r": sizes[b_type][g_type][0],  # Number of rows
                "c": sizes[b_type][g_type][1],  # Number of columns
            }

    # Generate the kernels
    kernel["K"] = {}
    for b_type in urban_challenges:
        kernel["K"][b_type] = {}
        for g_type in green_types:
            min_val, max_val = ranges[b_type][g_type]  # Get range for current green type and urban challenge
            rows, cols = sizes[b_type][g_type]  # Get size for current green type
            k = generate_kernel_2d(rows, cols, alg, min_val, max_val)
            kernel["K"][b_type][g_type] = np.array(k).flatten().tolist()  # Flatten kernel to a list

    return kernel


if __name__ == "__main__":
    # GENERATE KERNEL

    green_types = ["GreenWall", "GreenRoof", "StreetTree", "UrbanPark"]
    urban_challenges = ["TempMax", "TempMin", "Pm10", "Pm2", "Fairness"]

    alg = 3  # max to min

    ranges = {"TempMax": {"GreenWall": [0.1, 2.7], "GreenRoof": [0.1, 2.0], "StreetTree": [0.1, 1.3],
                          "UrbanPark": [0.1, 3.5]},
              "TempMin": {"GreenWall": [0.1, 1.9], "GreenRoof": [0.1, 1.4], "StreetTree": [0.1, 0.7],
                          "UrbanPark": [0.1, 2.5]},
              "Pm10": {"GreenWall": [0.01, 12.90], "GreenRoof": [0.01, 6.45], "StreetTree": [0.01, 10.32],
                       "UrbanPark": [0.01, 12.90]},
              "Pm2": {"GreenWall": [0.01, 5.03], "GreenRoof": [0.01, 2.51], "StreetTree": [0.01, 4.02],
                      "UrbanPark": [0.01, 5.03]},
              "Fairness": {"GreenWall": [2.0, 6.0], "GreenRoof": [0, 2.0], "StreetTree": [0.1, 4.0],
                           "UrbanPark": [4.0, 10.0]}
              }

    sizes = {"TempMax": {"GreenWall": [5, 5], "GreenRoof": [5, 5], "StreetTree": [5, 5],
                         "UrbanPark": [5, 5]},
             "TempMin": {"GreenWall": [3, 3], "GreenRoof": [3, 3], "StreetTree": [3, 3],
                         "UrbanPark": [3, 3]},
             "Pm10": {"GreenWall": [5, 5], "GreenRoof": [5, 5], "StreetTree": [3, 3],
                      "UrbanPark": [7, 7]},
             "Pm2": {"GreenWall": [5, 5], "GreenRoof": [5, 5], "StreetTree": [3, 3],
                     "UrbanPark": [7, 7]},
             "Fairness": {"GreenWall": [5, 5], "GreenRoof": [3, 3], "StreetTree": [3, 3],
                          "UrbanPark": [11, 11]}
             }

    kernel = generate_kernels_2d(alg, green_types, urban_challenges, sizes, ranges)

    with open("Kernels.json", 'w') as f:
        json.dump(kernel, f)
