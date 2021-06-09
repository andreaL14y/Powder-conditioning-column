# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
# import torchvision
# import torchvision.transforms as transforms

from torch.utils.data import DataLoader, TensorDataset

loss_function = nn.MSELoss()


class TooSimpleNetwork(nn.Module):

  def __init__(self, n_inputs, n_outputs):
    super().__init__()

    output0 = 8
    output1 = 4
    output2 = 2

    self.layer0 = nn.Linear(n_inputs, output0)
    #self.layer0 = nn.Linear(n_inputs, n_outputs)
    #self.layer1 = nn.Linear(output0, n_outputs)
    self.layer1 = nn.Linear(output0, output1)
    self.layer2 = nn.Linear(output1, output2)
    self.layer3 = nn.Linear(output2, n_outputs)

  def forward(self, x):
    # print('x in foward is', x)
    # print('x is type', type(x))
    # print('x is size', x.size())
    # x = x.view(-1, 1)
    # print('x is viewed to', x)
    x = self.layer0(x)
    # print('Done with layer 0')
    x = F.relu(x)
    x = self.layer1(x)
    x = F.relu(x)
    x = self.layer2(x)
    x = F.relu(x)
    x = self.layer3(x)

    output = x
    return output


def train_one_epoch(model, device, train_loader, optimizer, loss_function, scheduler, epoch):
  # print('Im in train_one_epoch')
  model.train()  # tells the model that it's training: This can alter some internal behavior
  log_interval = 129

  # print('Im in train_one_epoch')
  for batch_idx, (x, y) in enumerate(train_loader):
    # print('x is', x)
    # print('y is', y)
    x = x.to(device)
    y = y.to(device)

    optimizer.zero_grad()

    # execute the forward pass
    # print('About to predict something')
    prediction = model.forward(x)
    # print('Yeei, I predicted something :)')

    # calculate the loss
    prediction = prediction.view(-1, 1)
    y = y.view(-1, 1)
    loss = loss_function(prediction, y)
    loss.backward()
    optimizer.step()

    # print some updates
    if (batch_idx + 1) % log_interval == 0:
      print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.2e}'.format(
        epoch, (1 + batch_idx) * len(x), len(train_loader.dataset),
               100. * (1 + batch_idx) / len(train_loader), loss.item()))


def evaluate(model, device, loss_fn, val_loader):
  # Set internal model state to evaluation mode
  model.eval()
  val_loss = 0.0
  correct = 0

  # use the torch scope 'no_grad' to not calculate gradients in evaluation mode
  with torch.no_grad():
    for x, y in val_loader:
      x = x.to(device)
      y = y.to(device)
      prediction = model.forward(x)
      # prediction = prediction.view(-1, 1)
      # y = y.view(-1, 1)
      loss_values = loss_function(prediction, y)
      val_loss += loss_values.item()  # sum up batch loss

  val_loss /= len(val_loader.dataset)

  # and print some updates
  print('\nValidation set: Average loss: {:.2e}\n'.format(val_loss))
  return val_loss