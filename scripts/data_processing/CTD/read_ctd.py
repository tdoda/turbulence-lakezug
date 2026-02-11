# -*- coding: utf-8 -*-
import os
import sys
import json
import random
import logging
import shutil
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\ctd'))
from ctd import CTD
from functions_ctd import create_file_list, copy_files, read_data, process_profiles
logger = logging.getLogger(__name__)

#%% Specify field campaign here:

date_campaign='20260113'

#%% Other parameters

ctd_data_folder='..\..\data\ctd'
input_folder=os.path.join(ctd_data_folder,date_campaign)
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[
    logging.FileHandler("debug.log"),
    logging.StreamHandler()
])
extensions = [".rsk", ".tob", ".cnv", ".meta"]

for output_folder in ["Level1", "Failed"]:
    if os.path.exists(os.path.join(input_folder, output_folder)):
        print("Folder {} already exists: delete it".format(output_folder))
        shutil.rmtree(os.path.join(input_folder, output_folder))
    os.makedirs(os.path.join(input_folder, output_folder))

logger.info("Reprocessing data")
files = create_file_list(os.path.join(input_folder, "Level0"))

metadata_required = []

#%% Read and export CTD data
for file in files:
    logger.info("Processing file {}".format(file["path"]))
    metadata = True
    try:
        profiles = read_data(file["path"], file["type"])
    except Exception as e:
        print(e)
        logger.warning("Failed to process {}".format(file["path"]))
        copy_files(file["path"], os.path.join(input_folder, "Failed"), extensions)
        continue
    for profile in profiles:
        if os.path.isfile(os.path.join(os.path.dirname(file["path"]), profile["name"] + ".meta")):
            logger.info("Processing profile {}".format(profile["name"]))
            ctd = CTD(logger)
            if ctd.read_profile(profile):
                ctd.quality_assurance(r'..\..\functions\ctd\quality_assurance_ctd.json')
                lake = ctd.get_lake()
                file_name = os.path.basename(file["path"]).rsplit('.', 1)[0]
                ctd.export(os.path.join(input_folder, "Level1"), "L1_CTD_{}_{}".format(file["type"], file_name),overwrite=True)
                ctd.mask_data() # Replace flagged data by nan 
                ctd.derive_variables() # Compute additional variables to add to Level 2
                ctd.export(os.path.join(input_folder, "Level2"), "L2_CTD_{}_{}".format(file["type"], file_name), overwrite=True) # Create Level 2 file      
                
        else:
            logger.info("No metadata for profile {}".format(profile["name"]))
            metadata_required.append(profile)
            metadata = False


