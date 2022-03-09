#!/usr/bin/env python
# coding: utf-8



#libraries
import sys
import logging

import torch
import torch.nn as nn
import pandas as pd
import pickle
from collections import Counter
import numpy as np
from torchtext.data import Field 
from torchtext.data import Dataset, Example


# Libraries

import matplotlib.pyplot as plt
import pandas as pd
import torch
import torchtext
# Preliminaries

from torchtext.data import Field, TabularDataset, BucketIterator

# Models

import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence

# Training

import torch.optim as optim

# Evaluation

from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import seaborn as sns

logging.basicConfig(filename='app.log', filemode='w', format='%(asctime)s,%(msecs)d %(name)s - %(levelname)s - %(message)s',datefmt='%H:%M:%S', level=logging.DEBUG)

BATCH_SIZE = 128
MAX_VOCAB_SIZE = 160000
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def load_data():
    with open('./data/word_df_metrics.pkl', 'rb') as f:
        data = pickle.load(f)



def load_fields():
    # Fields for encoding
    tokenize = lambda x: x.split(' ')

    text_field = Field(sequential=True, tokenize=tokenize,lower=False, include_lengths=True, batch_first=True,pad_token='O')
    label_field = Field(sequential=False, use_vocab=False, batch_first=True, dtype=torch.float)

    #same order and same column name as in the dataframe 
    fields = [ ('label', label_field),('sequence_splitted', text_field)]

    return text_field, label_field, fields

def load_tabular_dataset(fields):    
    #train='train.csv', validation='valid.csv', test='test.csv',
    train, valid, test = TabularDataset.splits(path='./data/', train='valid1.csv', validation='test1.csv', test='test1.csv',
                                                format='CSV', fields=fields, skip_header=True)
    return train, valid, test

    # Iterators
def load_iterators(train, valid, test):   
    
    #each batch is of type torch.LongTensor, they are the numericalized batch
    train_iter = BucketIterator(train, batch_size=128, sort_key=lambda x: len(x.sequence_splitted),
                                device=DEVICE, sort=True, sort_within_batch=True)
    valid_iter = BucketIterator(valid, batch_size=128, sort_key=lambda x: len(x.sequence_splitted),
                                device=DEVICE, sort=True, sort_within_batch=True)
    test_iter = BucketIterator(test, batch_size=128, sort_key=lambda x: len(x.sequence_splitted),
                                device=DEVICE, sort=True, sort_within_batch=True)

    return train_iter, valid_iter, test_iter



class LSTM(nn.Module):

    def __init__(self, dimension, vocab_size, embedding_dim, hidden_dim, output_dim, n_layers, bidirectional, dropout, pad_idx):
        super().__init__()

        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx = pad_idx)
        self.dimension = dimension
        self.lstm = nn.LSTM(embedding_dim,
                            hidden_size=hidden_dim,
                            num_layers=n_layers,
                            batch_first=True,
                            bidirectional=bidirectional)
        self.drop = nn.Dropout(p=dropout)

        self.fc = nn.Linear(2*dimension, output_dim)

    def forward(self, text, text_len):

        text_emb = self.embedding(text)

        packed_input = pack_padded_sequence(text_emb, text_len, batch_first=True, enforce_sorted=False)
        packed_output, _ = self.lstm(packed_input)
        output, _ = pad_packed_sequence(packed_output, batch_first=True)

        out_forward = output[range(len(output)), text_len - 1, :self.dimension]
        out_reverse = output[:, 0, self.dimension:]
        out_reduced = torch.cat((out_forward, out_reverse), 1)
        text_fea = self.drop(out_reduced)

        text_fea = self.fc(text_fea)
        text_fea = torch.squeeze(text_fea, 1)
        text_out = torch.sigmoid(text_fea)

        
        return text_out


# Save and Load Functions

def save_checkpoint(save_path, model, optimizer, valid_loss):

    if save_path == None:
        return
    
    state_dict = {'model_state_dict': model.state_dict(),
                  'optimizer_state_dict': optimizer.state_dict(),
                  'valid_loss': valid_loss}
    
    torch.save(state_dict, save_path)
    logging.info(f'Model saved to ==> {save_path}')


def load_checkpoint(load_path, model, optimizer):

    if load_path==None:
        return
    
    state_dict = torch.load(load_path, map_location=DEVICE)
    logging.info(f'Model loaded from <== {load_path}')
    
    model.load_state_dict(state_dict['model_state_dict'])
    optimizer.load_state_dict(state_dict['optimizer_state_dict'])
    
    return state_dict['valid_loss']


def save_metrics(save_path, train_loss_list, valid_loss_list, global_steps_list):

    if save_path == None:
        return
    
    state_dict = {'train_loss_list': train_loss_list,
                  'valid_loss_list': valid_loss_list,
                  'global_steps_list': global_steps_list}
    
    torch.save(state_dict, save_path)
    logging.info(f'Model saved to ==> {save_path}')


def load_metrics(load_path):

    if load_path==None:
        return
    
    state_dict = torch.load(load_path, map_location=DEVICE)
    logging.info(f'Model loaded from <== {load_path}')
    
    return state_dict['train_loss_list'], state_dict['valid_loss_list'], state_dict['global_steps_list']




def train(model,optimizer,criterion ,train_loader):
    running_loss = 0.0 
    global_step = 0  
    model.train()   
    logging.info('train function')
    for (labels,(text, text_len)), _  in train_loader:           
            labels = labels.to(DEVICE)
            text = text.to(DEVICE)
            text_len = text_len.to(DEVICE)
            output = model(text, text_len)
 

            loss = criterion(output, labels)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # update running values
            running_loss += loss.item()
            # evaluate it only for few steps
            global_step += 1 
            if global_step % 20 == 0: 
                print(global_step)
    return running_loss, global_step  


def validate(model,criterion ,valid_loader):
    valid_running_loss = 0.0
    model.eval()
    logging.info('validate')
    with torch.no_grad():                    
        # validation loop
        for (labels, (text, text_len)), _ in valid_loader:
                    labels = labels.to(DEVICE)
                    text = text.to(DEVICE)
                    text_len = text_len.to(DEVICE)
                    output = model(text, text_len)

                    loss = criterion(output, labels)
                    valid_running_loss += loss.item()

    return valid_running_loss




def train_model(model, optimizer, criterion , train_loader , valid_loader , num_epochs , eval_every ,file_path ,
          best_valid_loss):
    print('start training')
    # initialize running values
    logging.info('train_model')
    global_step = 0
    train_loss_list = []
    valid_loss_list = []
    global_steps_list = []
    for epoch in range(num_epochs):
        print(epoch)
        train_loss,global_step = train(model,  optimizer, criterion,train_loader)
        logging.info('Finished Training!')
          #if global_step % eval_every == 0: I want to evaluate at each step
          #with torch.no_grad():
        val_loss = validate(model, criterion, valid_loader)
        logging.info('Finished validation!')  
        # evaluation
        average_train_loss = train_loss / len(train_loader) #eval_every
        average_valid_loss = val_loss / len(valid_loader)
        train_loss_list.append(average_train_loss)
        valid_loss_list.append(average_valid_loss)
        global_steps_list.append(global_step)
    
          # resetting running values
        train_loss = 0.0                
        val_loss = 0.0
        model.train()

        # print progress
        logging.info('Epoch [{}/{}], Step [{}/{}], Average Train Loss: {:.4f}, Average Valid Loss: {:.4f}, Train Loss: {:.4f}, Valid Loss: {:.4f}'
        .format(epoch+1, num_epochs, global_step, num_epochs*len(train_loader),
                              average_train_loss, average_valid_loss, train_loss, val_loss))
                
         # checkpoint
        if best_valid_loss > average_valid_loss:
            best_valid_loss = average_valid_loss
            logging.info(f'Best validation loss!! {best_valid_loss}')
            save_checkpoint(file_path + '/model.pt', model, optimizer, best_valid_loss)
            save_metrics(file_path + '/metrics.pt', train_loss_list, valid_loss_list, global_steps_list)
    
        save_metrics(file_path + '/metrics.pt', train_loss_list, valid_loss_list, global_steps_list)
        logging.info('Finished Training!')




def main():
    # print version and environment information
    print(f"{torch.__version__=}")
    print(f"{torch.version.cuda=}")
    print(f"{torch.backends.cudnn.enabled=}")
    print(f"{torch.cuda.is_available()=}")
    print(f"{DEVICE=}")
    if torch.cuda.is_available():
        print(f"{torch.cuda.get_device_properties(DEVICE)}")
       
    
    destination_folder='data/'
    
    text_field, label_field, fields = load_fields()

    train, valid, test=load_tabular_dataset(fields)

    train_iter, valid_iter, test_iter= load_iterators(train, valid, test)

    # Vocabulary
    text_field.build_vocab(train)
    label_field.build_vocab(train)

    logging.info(label_field.vocab.freqs)

    #Define model
    model = LSTM(dimension=128, vocab_size=len(text_field.vocab), embedding_dim=300, hidden_dim=128, output_dim=1, n_layers=1, 
                 bidirectional=True, dropout=0.5, pad_idx=text_field.vocab.stoi[text_field.pad_token]).to(DEVICE)


    optimizer = optim.Adam(model.parameters(), lr=0.001)
    train_model(model,optimizer,criterion = nn.BCELoss(),train_loader = train_iter,valid_loader = valid_iter,
          num_epochs = 1,eval_every = len(train_iter) // 2,file_path = destination_folder,
          best_valid_loss = float("Inf"))



if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print()
        print("Interrupted with CTRL-C, exiting...")
        sys.exit()

