
# get the index of the machine generated images using the precision of the angle 
# [ https://www.kaggle.com/brassmonkey381/viewing-leak-and-machine-images ,
#   https://www.kaggle.com/c/statoil-iceberg-classifier-challenge/discussion/46235#261557 ]


import numpy as np
import pandas as pd
import re

PATH_TEST = "../test.json"

dat_test = pd.read_json(PATH_TEST, precise_float = True)

angle_test = dat_test.iloc[:, 3]

angle_test = np.array(angle_test.as_matrix(columns=None))

#-----------------------------------------------------------------------------------------------------------------
# first convert angle to string then split the string and count the decimal points of the second part of the list
# if an angle item has 4 decimal places then it is a true image ( NOT a machine generated image )
#-----------------------------------------------------------------------------------------------------------------

angle_prec = [str(angle_test[i]) for i in range(len(angle_test))]

angle_spl = [re.split('[.]', str(angle_prec[i])) for i in range(len(angle_prec))]

angle_count = [len(angle_spl[i][1]) for i in range(len(angle_spl))]

bool_true_images = [i <= 4 for i in angle_count]                                 # some numbers have less than 4 after period digits

print(sum(bool_true_images))                                                     # 3424 images in total

idx_true_images = [i for i in range(len(bool_true_images)) if bool_true_images[i] == True]

pd_df = pd.DataFrame()

pd_df['true_images'] = idx_true_images

pd_df.to_csv("../true_indices_test_data.csv", index=False)
