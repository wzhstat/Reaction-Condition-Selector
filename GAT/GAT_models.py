import argparse
import os.path as osp

import torch
import torch.nn.functional as F

from torch.nn import Linear, Sequential, Dropout
import torch_geometric.transforms as T
from torch_geometric.datasets import Planetoid
from torch_geometric.logging import init_wandb, log
from torch_geometric.nn import GATConv
from torch_geometric.nn import global_mean_pool

class GAT_test(torch.nn.Module):
    def __init__(self, num_features, num_classes, num_layers=3, num_heads=8, dropout=0.2, edge_channels = 147, add_out_ffn = False):
        super(GAT_test, self).__init__()
        self.dropout = dropout
        self.add_out_ffn = add_out_ffn

        # Create multiple layers of GAT
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers):
            in_channels = num_features if i == 0 else num_heads * num_features
            out_channels = num_features if i < num_layers - 1 else num_classes
            self.convs.append(GATConv(in_channels, out_channels, heads=num_heads, dropout=dropout, edge_dim = edge_channels))
        self.classifier = Linear(num_heads*out_channels, num_classes)
        self.edge_ffn = Sequential(
            Linear(edge_channels, edge_channels),
            Dropout(dropout),
        )
        self.atom_ffn = Sequential(
            Linear(num_features, num_features),
            Dropout(dropout),
        )
        self.out_ffn = Sequential(
            Linear(num_classes,num_classes),
            Dropout(dropout),
        )
        

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch

        if self.add_out_ffn:
            x = self.atom_ffn(x)
            edge_attr = self.edge_ffn(edge_attr)

        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index,edge_attr)
            if i < len(self.convs) - 1:
                x = F.elu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        
        unique_batches = batch.unique(sorted=True)
        batch_mapped = torch.zeros_like(batch)
        for i, unique_batch in enumerate(unique_batches):
            batch_mapped[batch == unique_batch] = i
        x = global_mean_pool(x, batch_mapped)  # Global mean pooling
        x = self.classifier(x)
        if self.add_out_ffn:
            x = self.out_ffn(x)
        return x

class GAT(torch.nn.Module):
    def __init__(self, num_features, num_classes, num_layers=3, num_heads=8, dropout=0.2, edge_channels = 147):
        super(GAT, self).__init__()
        self.dropout = dropout

        # Create multiple layers of GAT
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers):
            in_channels = num_features if i == 0 else num_heads * num_features
            out_channels = num_features if i < num_layers - 1 else num_classes
            self.convs.append(GATConv(in_channels, out_channels, heads=num_heads, dropout=dropout, edge_dim = edge_channels))
        self.classifier = Linear(num_heads*out_channels, num_classes)

        

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        '''
        if self.add_in_ffn:
            x = self.atom_ffn(x)
            edge_attr = self.edge_ffn(edge_attr)
        '''

        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index,edge_attr)
            if i < len(self.convs) - 1:
                x = F.elu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        
        unique_batches = batch.unique(sorted=True)
        batch_mapped = torch.zeros_like(batch)
        for i, unique_batch in enumerate(unique_batches):
            batch_mapped[batch == unique_batch] = i
        x = global_mean_pool(x, batch_mapped)  # Global mean pooling
        x = self.classifier(x)
        return x