import  numpy as np

import  torch
from    torch import nn
from    torch import optim
from    torch.nn import functional as F


from    data import load_data, preprocess_features, preprocess_adj, sample_mask
import  model
from    config import  args
from    utils import masked_loss, masked_acc
import  pickle as pkl
import  scipy.sparse as sp

import random

seed = 123
np.random.seed(seed)
torch.random.manual_seed(seed)


if torch.cuda.is_available():
    torch.cuda.set_device(0)
else:
    print("Running with cpu")


adj        = pkl.load(open("Cyber_data/contig.graph",'rb'))
labels     = pkl.load(open("Cyber_data/contig.label",'rb'))
features   = pkl.load(open("Cyber_data/contig.feature",'rb'))
test_to_id = pkl.load(open("Cyber_data/contig.dict",'rb'))
idx_test   = pkl.load(open("Cyber_data/contig.mask",'rb'))

idx_test = np.array(idx_test)
labels = np.array(labels)

y_train = np.zeros(labels.shape)
y_test = np.zeros(labels.shape)



idx_train = np.array([i for i in range(len(labels)) if i not in idx_test])


train_mask = sample_mask(idx_train, labels.shape[0])
test_mask = sample_mask(idx_test, labels.shape[0])

y_train[train_mask] = labels[train_mask]
y_test[test_mask] = labels[test_mask]


features = sp.csc_matrix(features)

print('adj:', adj.shape)
print('features:', features.shape)
print('y:', y_train.shape, y_test.shape) # y_val.shape, 
print('mask:', train_mask.shape, test_mask.shape) # val_mask.shape

features = preprocess_features(features) # [49216, 2], [49216], [2708, 1433]
supports = preprocess_adj(adj)

if torch.cuda.is_available():
    torch.cuda.set_device(0)
    device = torch.device('cuda')
    train_label = torch.from_numpy(y_train).long().to(device)
    num_classes = max(labels)+1
    train_mask = torch.from_numpy(train_mask.astype(np.bool)).to(device)
    test_label = torch.from_numpy(y_test).long().to(device)
    test_mask = torch.from_numpy(test_mask.astype(np.bool)).to(device)

    i = torch.from_numpy(features[0]).long().to(device)
    v = torch.from_numpy(features[1]).to(device)
    feature = torch.sparse.FloatTensor(i.t(), v, features[2]).float().to(device)

    i = torch.from_numpy(supports[0]).long().to(device)
    v = torch.from_numpy(supports[1]).to(device)
    support = torch.sparse.FloatTensor(i.t(), v, supports[2]).float().to(device)

else:
    train_label = torch.from_numpy(y_train).long()
    num_classes = max(labels)+1
    train_mask = torch.from_numpy(train_mask.astype(np.bool))
    test_label = torch.from_numpy(y_test).long()
    test_mask = torch.from_numpy(test_mask.astype(np.bool))

    i = torch.from_numpy(features[0]).long()
    v = torch.from_numpy(features[1])
    feature = torch.sparse.FloatTensor(i.t(), v, features[2]).float()

    i = torch.from_numpy(supports[0]).long()
    v = torch.from_numpy(supports[1])
    support = torch.sparse.FloatTensor(i.t(), v, supports[2]).float()

print('x :', feature)
print('sp:', support)
num_features_nonzero = feature._nnz()
feat_dim = feature.shape[1]


def accuracy(out, mask):
    pred = np.argmax(out, axis = 1)
    mask_pred = np.array([pred[i] for i in range(len(labels)) if mask[i] == True])
    mask_label = np.array([labels[i] for i in range(len(labels)) if mask[i] == True])
    return np.sum(mask_label == mask_pred)/len(mask_pred)

net = model.GCN(feat_dim, num_classes, num_features_nonzero)
if torch.cuda.is_available():
    net.to(device)
optimizer = optim.Adam(net.parameters(), lr=0.01)#args.learning_rate


_ = net.train()
for epoch in range(args.epochs*2):
    # forward pass
    out = net((feature, support))
    #out = out[0]
    loss = masked_loss(out, train_label, train_mask)
    loss += args.weight_decay * net.l2_loss()
    # backward pass
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    # output
    if epoch % 10 == 0:
        # calculating the acc
        _ = net.eval()
        out = net((feature, support))
        if torch.cuda.is_available():
            acc_train = accuracy(out.detach().cpu().numpy(), train_mask.detach().cpu().numpy())
        else:
            acc_train = accuracy(out.detach().numpy(), train_mask.detach().numpy())
        #acc_test = accuracy(out.detach().cpu().numpy(), test_mask.detach().cpu().numpy())
        print(epoch, loss.item(), acc_train)
        if acc_train > 0.978:
            break
    _ = net.train()


net.eval()
out = net((feature, support))
out = F.softmax(out,dim =1)
if torch.cuda.is_available():
    out = out.cpu().detach().numpy()
else:
    out = out.detach().numpy()

pred = np.argmax(out, axis = 1)
score = np.max(out, axis = 1)

mode = "testing"
if mode == "validation":
    print(classification_report(labels, pred))
    print(accuracy(out, train_mask.detach().cpu().numpy()))
    mask = test_mask.detach().cpu().numpy()
    test_pred = np.array([pred[i] for i in range(len(pred)) if mask[i] == True])
    test_label = np.array([labels[i] for i in range(len(labels)) if mask[i] == True])
    print(classification_report(test_label, test_pred))
    print(np.sum(test_label == test_pred)/len(test_pred))


pred_to_label = {0:"Ackermannviridae", 1:"Autographiviridae", 2:"Demerecviridae",
3:"Drexlerviridae", 4:"Herelleviridae", 5:"Myoviridae", 6:"Podoviridae", 7:"Siphoviridae"}


with open("prediction.csv", 'w') as f_out:
    _ = f_out.write("contig_names,prediction,score\n")
    for key in test_to_id.keys():
        if labels[test_to_id[key]] == -1:
            _ = f_out.write(str(key) + "," + str(pred_to_label[pred[test_to_id[key]]]) + "," + str(score[test_to_id[key]]) + "\n")
        else:
            _ = f_out.write(str(key) + "," + str(pred_to_label[labels[test_to_id[key]]]) + "," + str(1) + "\n")

