import torch
import numpy as np
from torch import nn
from rdkit import Chem

# define highway layer
class highway(nn.Module):
    def __init__(self,size):
        super().__init__()
        self.H = nn.Linear(size, size)
        self.T = nn.Linear(size, size)
        self.T.bias.data.fill_(-1.)
    def forward(self,x):
        h = torch.relu(self.H(x))
        t = torch.sigmoid(self.T(x))
        c = 1. - t
        y = h * t + x * c
        return y

class nnModel0(nn.Module):
    '''
    A MLP model for predicting the reaction conditions.
    input: 512-dim reaction fingerprint, 512-dim product fingerprint
    output: x-dim reaction conditions
    The input pass through two nnlayers.
    '''
    def __init__(self, out_size, n0 ,n1 = 128, n2 = 32 ,dp = 0.1):
        super(nnModel0,self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(n0,n1),
            nn.Dropout(dp),
            nn.ELU()
        )
        self.l2 = nn.Sequential(
            nn.Linear(n1,n2),
            nn.Dropout(dp),
            nn.ELU()
        )
        self.out = nn.Sequential(
            nn.Linear(n2,out_size),
            nn.Softmax()
        )
    def forward(self,xs):
        x = torch.cat(xs,dim=1)
        x = self.l2(self.l1(x))
        return self.out(x)


class nnModel1(nn.Module):
    '''
    A MLP model for predicting the reaction conditions.
    input: 512-dim reaction fingerprint, 512-dim product fingerprint
    output: x-dim reaction conditions
    The input pass through a nn layer and two highway layer.
    '''
    def __init__(self, out_size,n0 ,n1 = 128, n2 = 32 ,dp = 0.1):
        super(nnModel1,self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(n0,n1),
            nn.Dropout(dp),
            nn.ELU()
        )

        self.l2 = [nn.Sequential(  
            highway(n1),
            nn.ELU()
        ) for _ in range(2)]

        self.out = nn.Sequential(
            nn.Linear(n1,out_size),
            nn.Softmax()
        )
    def forward(self,xs):
        x = torch.cat(xs,dim=1)
        x3 = self.l1(x)
        for layer in self.l2:
            x3 = layer(x3)
        return self.out(x3)


class nnModel2(nn.Module):
    '''
    A MLP model for predicting the reaction conditions.
    input: 512-dim reaction fingerprint, 512-dim product fingerprint, n-dim one-hot template
    output: x-dim reaction conditions
    The input pass through a nn layer and two highway layer.
    '''
    def __init__(self, out_size, n, n1=128, n2=32, dp=0.1):
        super(nnModel2, self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(1024 + n, n1),
            nn.Dropout(dp),
            nn.ELU()
        )

        self.l2 = [nn.Sequential(
            highway(n1),
            nn.ELU()
        ) for _ in range(2)]

        self.out = nn.Sequential(
            nn.Linear(n1, out_size),
            nn.Dropout(dp),
            nn.Softmax()
        )

    def forward(self, xs):
        x = torch.cat(xs, dim=1)
        x = self.l1(x)
        for layer in self.l2:
            x = layer(x)
        x = self.out(x)
        return x




