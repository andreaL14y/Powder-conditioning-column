# IMPORTS
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
# import torchvision
# import torchvision.transforms as transforms

from torch.utils.data import DataLoader, TensorDataset

from nn_functions import TooSimpleNetwork, evaluate, train_one_epoch, loss_function

######################################### LOAD ALL DATA ################################################################
am_material_40_24         = np.load('RH-0.4_T-24.0_AmMaterial.npy')
print(am_material_40_24.shape)
ttg_material_40_24        = np.load('RH-0.4_T-24.0_TTg.npy')
am_material_40_24_new     = am_material_40_24[1:, :, :]
print(am_material_40_24_new.shape)
am_material_40_24         = am_material_40_24[:-1, :, :]
ttg_material_40_24        = ttg_material_40_24[:-1, :, :]

am_material_50_35         = np.load('RH-0.5_T-35.0_AmMaterial.npy')
ttg_material_50_35        = np.load('RH-0.5_T-35.0_TTg.npy')
am_material_50_35_new     = am_material_50_35 [1:, :, :] #- am_material_50_35[:-1, :, :]
am_material_50_35         = am_material_50_35 [:-1, :, :]
ttg_material_50_35        = ttg_material_50_35[:-1, :, :]

am_material_60_24         = np.load('RH-0.6_T-297.15_AmMaterial.npy')
ttg_material_60_24        = np.load('RH-0.6_T-297.15_TTg.npy')
am_material_60_24_new     = am_material_60_24 [1:, :, :] #- am_material_60_24[:-1, :, :]
am_material_60_24         = am_material_60_24 [:-1, :, :]
ttg_material_60_24        = ttg_material_60_24[:-1, :, :]

am_material_70_35         = np.load('RH-0.7_T-35.0_AmMaterial.npy')
ttg_material_70_35        = np.load('RH-0.7_T-35.0_TTg.npy')
am_material_70_35_new     = am_material_70_35[1:, :, :] #- am_material_70_35[:-1, :, :]
am_material_70_35         = am_material_70_35[:-1, :, :]
ttg_material_70_35        = ttg_material_70_35[:-1, :, :]

am_material_40_24         = am_material_40_24.flatten()
ttg_material_40_24        = ttg_material_40_24.flatten()
am_material_40_24_new     = am_material_40_24_new.flatten()

am_material_50_35         = am_material_50_35.flatten()
ttg_material_50_35        = ttg_material_50_35.flatten()
am_material_50_35_new     = am_material_50_35_new.flatten()

am_material_60_24         = am_material_60_24.flatten()
ttg_material_60_24        = ttg_material_60_24.flatten()
am_material_60_24_new     = am_material_60_24_new.flatten()

am_material_70_35         = am_material_70_35.flatten()
ttg_material_70_35        = ttg_material_70_35.flatten()
am_material_70_35_new     = am_material_70_35_new.flatten()


################################# CHOOSE AND PROCESS DATA ##############################################################
am_material = am_material_50_35
ttg_material = ttg_material_50_35
am_material_new = am_material_50_35_new

# WANT TO INPUT AMOUNT AM MATERIAL AND T-Tg, GET AM MATERIAL t+1
x = torch.tensor([am_material, ttg_material], dtype=torch.float32)
x = torch.transpose(x, 0, 1)

y = torch.tensor(am_material_new, dtype=torch.float32)
print('x has shape:', x.shape)
print('y has shape:', y.shape)
# start = 9005
# print(np.around(x[start:start+15], 3))
# print(np.around(y[start:start+15], 3))

# Set proportions of data, and compute number of data points of each feature
n_data_points = len(x)
train = 0.55
val = 0.25
test = 1 - train - val

n_train = int(round(n_data_points * train))
n_val = int(round(n_data_points * val))
n_test = int(round(n_data_points - n_train - n_val))
if n_test + n_train + n_val != n_data_points:
  print('Some numbers are messed up :(')

# Create dataset of all info and split on randomized indeces
full_dataset = TensorDataset(x, y)
indeces = np.arange(0, n_data_points)
[train_indeces, val_indeces, test_indeces] = torch.utils.data.random_split(indeces, [n_train, n_val, n_test])

#train_indeces = list(train_indeces)

# Create separate datasets and get the correct dimensions
train_dataset = full_dataset[train_indeces]
train_dataset = list(zip(*train_dataset))

val_dataset = full_dataset[val_indeces]
val_dataset = list(zip(*val_dataset))

test_dataset = full_dataset[test_indeces]
test_dataset = list(zip(*test_dataset))

print(f'There are {n_data_points} data points in total')
print(train_dataset[0:2])

################################## Set up parameters for training ######################################################
n_inputs = 2                      # current amount amorph and T-Tg
n_outputs = 1                     # prediction of next amount amorph

model = TooSimpleNetwork(n_inputs, n_outputs)

# Use a GPU if it is available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
cpu_device = torch.device("cpu")
print('Using', device, 'to compute')
model = model.to(device)

# Only with this information can pytorch work its magic
alpha = 0.0015
optimizer = optim.Adam(model.parameters(), lr = alpha)

# And optionally, a learning rate scheduler
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer)

####################################### TRAINING #######################################################################
from copy import deepcopy

batch_size = 32

print('Batch size is:', batch_size)
print(f'There will be {len(train_dataset)/batch_size:.1f} batches in total')

train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
print('Train loader is', len(train_loader), 'long\n')
val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size)

# 20 epochs
num_epochs = 20

best_model = None
best_val_loss = 9999999

for epoch in range(1, num_epochs + 1):
    # Training and Evaluation
    train_one_epoch(model, device, train_loader, optimizer, loss_function, scheduler, epoch)
    val_loss = evaluate(model, device, loss_function, val_loader)

    scheduler.step(val_loss)

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        best_model = deepcopy(model.to(cpu_device))
        model.to(device)

########################################### TEST #######################################################################
model.eval()

test_dataset_list = list(zip(*test_dataset))
test_dataset_list_x = list(test_dataset_list[0])
test_dataset_list_y = list(test_dataset_list[1])
test_dataset_tensor = torch.stack(test_dataset_list_x)

with torch.no_grad():
  test_predictions = model(test_dataset_tensor)

# start_n = 7000
# print(np.around(test_dataset_tensor[start_n:start_n+15], 3))
# print(np.around(test_predictions[start_n:start_n+15], 3))

test_predictions = test_predictions.numpy()
test_predictions = test_predictions.squeeze()

plt.plot(test_predictions, test_dataset_list_y, '+', color='orange', label='Predictions')
plt.xlabel('Predicted by nn')
plt.ylabel('True value')
plt.legend()
plt.show()