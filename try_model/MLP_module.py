import torch
import numpy as np
from torch import nn
from rdkit import Chem

class MLP_module0(nn.Module):
    def __init__(self,out_size):
        super(MLP_module0,self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(512,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l2 = nn.Sequential(
            nn.Linear(128,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l3 = nn.Sequential(
            nn.Linear(512,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l4 = nn.Sequential(
            nn.Linear(128,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l5 = nn.Sequential(
            nn.Linear(16,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l6 = nn.Sequential(
            nn.Linear(128,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l7 = nn.Sequential(
            nn.Linear(128,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.out = nn.Sequential(
            nn.Linear(32,out_size),
            nn.Softmax()
        )

    def forward(self,rfpgen,pfpgen):
        x1 = self.l2(self.l6(self.l1(rfpgen)))
        x2 = self.l4(self.l6(self.l3(pfpgen)))
        return self.out(x1*x2)

class MLP_module1(nn.Module):
    def __init__(self,out_size):
        super(MLP_module1,self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(512,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l2 = nn.Sequential(
            nn.Linear(128,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l3 = nn.Sequential(
            nn.Linear(512,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l4 = nn.Sequential(
            nn.Linear(128,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l6 = nn.Sequential(
            nn.Linear(128,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l7 = nn.Sequential(
            nn.Linear(128,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.out = nn.Sequential(
            nn.Linear(32,out_size),
            nn.Softmax()
        )

    def forward(self,rfpgen,rxnfp):
        x1 = self.l2(self.l6(self.l1(rfpgen)))
        x2 = self.l4(self.l6(self.l3(rxnfp)))
        return self.out(x1*x2)




class MLP_module2(nn.Module):
    def __init__(self,out_size):
        super(MLP_module2,self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(512,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l2 = nn.Sequential(
            nn.Linear(128,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l3 = nn.Sequential(
            nn.Linear(512,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l4 = nn.Sequential(
            nn.Linear(128,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l5 = nn.Sequential(
            nn.Linear(16,32),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l6 = nn.Sequential(
            nn.Linear(128,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l7 = nn.Sequential(
            nn.Linear(128,128),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.out = nn.Sequential(
            nn.Linear(32,out_size),
            nn.Softmax()
        )

    def forward(self,rfpgen,pfpgen,b_tem):
        x1 = self.l2(self.l6(self.l1(rfpgen)))
        x2 = self.l4(self.l6(self.l3(pfpgen)))
        x3 =  self.l5(b_tem)
        return self.out((x1*x2)+x3*2)



class MLP_module3(nn.Module):
    def __init__(self,out_size):
        super(MLP_module3,self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(512,512),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l3 = nn.Sequential(
            nn.Linear(512,512),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l5 = nn.Sequential(
            nn.Linear(25065,512),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l6 = nn.Sequential(
            nn.Linear(512,512),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.l7 = nn.Sequential(
            nn.Linear(512,512),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.out = nn.Sequential(
            nn.Linear(512,out_size),
            nn.Softmax()
        )

    def forward(self,rfpgen,pfpgen,b_tem):
        x1 =self.l6(self.l1(rfpgen))
        x2 =self.l6(self.l3(pfpgen))
        x3 =self.l5(b_tem)
        return self.out((x1*x2)+x3*5)


class MLP_module5(nn.Module):
    def __init__(self,out_size):
        super(MLP_module5,self).__init__()
        self.l1 = nn.Sequential(
            nn.Linear(25065,512),
            nn.Dropout(0.2),
            nn.ELU()
        )
        self.out = nn.Sequential(
            nn.Linear(512,out_size),
            nn.Softmax()
        )
    def forward(self,b_tem):
        return  self.out(self.l1(b_tem))


