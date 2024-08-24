import torch
from GAT.GAT_models import GAT
import torch_geometric
from GAT.get_dataset import GraphDataset, collate_fn, GraphDataLoader
import pandas as pd
import os 


if torch.cuda.is_available():
    device = torch.device('cuda')
elif torch_geometric.is_xpu_available():
    device = torch.device('xpu')
else:
    device = torch.device('cpu')


def train(model, dataloader, val_dataloader, criterion, optimizer, num_epochs = 30, device = 'cuda', target = 'cat'):
    model.train()
    model.to(device)
    evals = eval(model, val_dataloader, device, criterion)
    last_lost = evals["Loss"]
    for epoch in range(num_epochs):
        total_loss = 0
        for i, batch in enumerate(dataloader):
            optimizer.zero_grad()
            batch = batch.to(device)
            logits = model(batch)
            loss = criterion(logits, batch.y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            if i % 1000 == 0:
                print(f'Epoch {epoch+1}, Loss: {total_loss / (i+1)}')
        print(f'Epoch {epoch+1}, Loss: {total_loss / len(dataloader)}')
        evals = eval(model, val_dataloader, device, criterion)
        if evals["Loss"] < last_lost:
            last_lost = evals["Loss"]
            if not os.path.exists('./GATmodels'):
                os.makedirs('./GATmodels')
            torch.save(model.state_dict(), './GATmodels/model_%s.pth'%(target))
    return model

def eval(model, dataloader, device, criterion):
    model.eval()
    model.to(device)
    with torch.no_grad():
        total_loss = 0
        top1_correct = 0
        top3_correct = 0
        top5_correct = 0
        top10_correct = 0
        total_samples = 0

        for batch in dataloader:
            batch = batch.to(device)
            outputs = model(batch)
            loss = criterion(outputs, batch.y)
            total_loss += loss.item()

            # Get the top 10 predictions from the logits
            _, pred = outputs.topk(10, 1, True, True)
            pred = pred.t()
            correct = pred.eq(batch.y.view(1, -1).expand_as(pred))

            # Calculate Top-1, Top-3, Top-5, and Top-10 accuracies
            top1_correct += correct[:1].reshape(-1).float().sum(0, keepdim=True).item()
            top3_correct += correct[:3].reshape(-1).float().sum(0, keepdim=True).item()
            top5_correct += correct[:5].reshape(-1).float().sum(0, keepdim=True).item()
            top10_correct += correct[:10].reshape(-1).float().sum(0, keepdim=True).item()

            total_samples += batch.y.size(0)

        print(f'Top-1 Accuracy: {top1_correct / total_samples:.4f}')
        print(f'Top-3 Accuracy: {top3_correct / total_samples:.4f}')
        print(f'Top-5 Accuracy: {top5_correct / total_samples:.4f}')
        print(f'Top-10 Accuracy: {top10_correct / total_samples:.4f}')
        print(f'Loss: {total_loss / len(dataloader):.4f}')

    return {
        'Top-1 Accuracy': top1_correct / total_samples,
        'Top-3 Accuracy': top3_correct / total_samples,
        'Top-5 Accuracy': top5_correct / total_samples,
        'Top-10 Accuracy': top10_correct / total_samples,
        'Loss': total_loss / len(dataloader)
    }


if __name__ == '__main__':  
    data_train = pd.read_csv('./GCN_data/GCN_data_train.csv')
    data_test = pd.read_csv('./GCN_data/GCN_data_test.csv')
    data_val = pd.read_csv('./GCN_data/GCN_data_val.csv')

    target_list = ['cat','solv0',	'solv1',	'reag0',	'reag1',	'reag2']
    target_class_num = [439, 542, 542, 2746, 2746, 2746]
    for i in range(6):
        print('Training for target: ', target_list[i])
        print('Class number: ', target_class_num[i])
        target = target_list[i]
        # get test dataset
        test_dataset = GraphDataset(data_test,target)
        test_loader = GraphDataLoader(test_dataset, batch_size=16, collate_fn=collate_fn)
        # get val dataset
        val_dataset = GraphDataset(data_val,target)
        val_loader = GraphDataLoader(val_dataset, batch_size=16, collate_fn=collate_fn)
        # get train dataset
        train_dataset = GraphDataset(data_train,target)
        train_loader = GraphDataLoader(train_dataset, batch_size=16, collate_fn=collate_fn)
        # Define the model
        node_in_channels = 165
        edge_in_channels = 193
        num_classes = target_class_num[i]
        hidden_channels = 1024
        num_layers = 5
        num_heads = 8
        dropout =0.2
        model = GAT(node_in_channels, num_classes, num_layers, num_heads, dropout, edge_in_channels)

        # Define the loss and optimizer
        criterion = torch.nn.CrossEntropyLoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=0.0001, weight_decay=1e-4, amsgrad = True)

        lambda1 = lambda epoch: 0.95 ** epoch
        scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer,lr_lambda = lambda1)

        # Train the model
        train(model, train_loader, val_loader, criterion, optimizer, num_epochs = 60, device='cuda', target=target)

        # Evaluate the model
        eval(model, test_loader, 'cuda', criterion)
