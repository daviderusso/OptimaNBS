green_types = ["GreenWall", "GreenRoof", "StreetTree", "UrbanPark"]  # NBS types
urban_challenges = ["TempMax", "TempMin", "Pm10", "Pm2", "Fairness"]  # Urban Challenge considered

alg_type = 3  # alg 3 is from max value in the center of the kernel to min values on the border

max_val_default = 1.0  # kernel default max value
min_val_default = 0.1  # kernel default min value

# ranges used to create kernels for each Urban Challenge and NBS type
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

# sizes of the kernels created for each Urban Challenge and NBS type
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

output_name = "Kernels.json"  # output filename
