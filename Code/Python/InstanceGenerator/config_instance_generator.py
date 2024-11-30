folderBase = '../Data/'
folderData = 'NationalData/'
folderMunicipality = 'Shp_Municipality/'
folderTiles = '/Tiles/'
folderInstance = '/Instances/'

image_format = '.tif'
shp_format = '.shp'
instance_format = ".json"

raw = "_raw"

cities = ['Bologna', 'Catania', 'Milan', 'Naples', 'Rome', 'Turin']
# cities = ['Bologna']

# TempMax and TempMin data refer to the 19/07/
dataName = ["LandUseGHS", "TempMax", "TempMin", "Pm10", "Pm2", "Fairness"]

# TempMax and TempMin data are multiplied by 100. So to report value to their normal magnitude
# are divided by 100 with the data normalizer
data_normalizer = dict(LandUseGHS=1, TempMax=0.01, TempMin=0.01, Pm10=1, Pm2=1,
                       Fairness=1)

# Reported in the paper as  XS, S, M, L, respectively
splitSize = [50, 100, 200, 300]
splitLabel = ["XS", "S", "M", "L"]

perc_to_be_valid = 0.25

green_type = ["GreenWall", "GreenRoof", "StreetTree", "UrbanPark"]

# LANDUSE CONFIG
forbidden_for_all = [4, 5, 21, 22, 23, 24, 25, 255]
allowed_roof = [14, 15]
allowed_wall = [11, 12, 13]
allowed_for_all = [1, 3]
pre_existent = {"GreenWall": [], "GreenRoof": [], "StreetTree": [2], "UrbanPark": []}

urban_challenge = ["TempMax", "TempMin", "Pm10", "Pm2", "Fairness"]

# INSTANCES CREATION JSON STRING
W = "W"
H = "H"
GT = "GreenType"
UC = "UrbanChallenge"
Forbid = "Forbidden"
PreExist = "PreExistent"
A = "A"

