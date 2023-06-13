import numpy as np
from Bio import SeqIO
from nltk import trigrams, bigrams
from gensim.models import Word2Vec
import re

np.set_printoptions(threshold=np.inf)


def w2v_fea_hm(fa_file, w2v_model):
    texts = []
    for index, record in enumerate(SeqIO.parse(fa_file, 'fasta')):
        tri_tokens = bigrams(record.seq)
        temp_str = ""
        for item in ((tri_tokens)):
            # print(item)
            temp_str = temp_str + " " + item[0] + item[1]
            # temp_str = temp_str + " " +item[0]
        texts.append(temp_str)

    seq = []
    stop = '[’!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+'
    for doc in texts:
        doc = re.sub(stop, '', doc)
        seq.append(doc.split())

    w2v_model = Word2Vec.load(w2v_model)

    embedding_matrix = w2v_model.wv.vectors

    ze = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    embedding_matrix = np.row_stack((embedding_matrix, ze))

    vocab_list = list(w2v_model.wv.key_to_index.keys())

    word_index = {word: index for index, word in enumerate(vocab_list)}

    def get_index(sentence):
        sequence = []

        ze_num = [426]
        for word in sentence:
            if word in vocab_list:
                sequence.append(word_index[word])
            else:
                sequence.extend(ze_num)
        return sequence

    X_data = np.array(list(map(get_index, seq)))  # shape (n, 28)

    vector = np.hstack((embedding_matrix[X_data][:,0,:], embedding_matrix[X_data][:,1,:], embedding_matrix[X_data][:,2,:]
                        , embedding_matrix[X_data][:,3,:], embedding_matrix[X_data][:,4,:], embedding_matrix[X_data][:,5,:]
                        , embedding_matrix[X_data][:,6,:], embedding_matrix[X_data][:,7,:], embedding_matrix[X_data][:,8,:]
                        , embedding_matrix[X_data][:,9,:], embedding_matrix[X_data][:,10,:], embedding_matrix[X_data][:,11,:]
                        , embedding_matrix[X_data][:,12,:], embedding_matrix[X_data][:,13,:], embedding_matrix[X_data][:,14,:]
                        , embedding_matrix[X_data][:,15,:], embedding_matrix[X_data][:,16,:], embedding_matrix[X_data][:,17,:]
                        , embedding_matrix[X_data][:,18,:], embedding_matrix[X_data][:,19,:], embedding_matrix[X_data][:,20,:]
                        , embedding_matrix[X_data][:,21,:], embedding_matrix[X_data][:,22,:], embedding_matrix[X_data][:,23,:]
                        , embedding_matrix[X_data][:,24,:], embedding_matrix[X_data][:,25,:], embedding_matrix[X_data][:,26,:]
                        , embedding_matrix[X_data][:,27,:]))  # shape (n, 840)

    return vector


def w2v_fea_ms(fa_file, w2v_model):
    texts = []
    for index, record in enumerate(SeqIO.parse(fa_file, 'fasta')):
        tri_tokens = bigrams(record.seq)
        temp_str = ""
        for item in ((tri_tokens)):
            # print(item)
            temp_str = temp_str + " " + item[0] + item[1]
            # temp_str = temp_str + " " +item[0]
        texts.append(temp_str)

    seq = []
    stop = '[’!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+'
    for doc in texts:
        doc = re.sub(stop, '', doc)
        seq.append(doc.split())

    # print(seq)  # 分割后的二元氨基酸  441个

    w2v_model = Word2Vec.load(w2v_model)

    embedding_matrix = w2v_model.wv.vectors

    ze = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    embedding_matrix = np.row_stack((embedding_matrix, ze))

    vocab_list = list(w2v_model.wv.key_to_index.keys())
    # print(vocab_list)   # voctor 分割的二元氨基酸 426个

    word_index = {word: index for index, word in enumerate(vocab_list)}
    # print(word_index)  # 分割的二元氨基酸索引 0-425

    def get_index(sentence):
        sequence = []

        ze_num = [425]
        for word in sentence:
            if word in vocab_list:
                sequence.append(word_index[word])
            else:
                sequence.extend(ze_num)
        return sequence


    X_data = np.array(list(map(get_index, seq)))  # shape (n, 28)

    vector = np.hstack((embedding_matrix[X_data][:,0,:], embedding_matrix[X_data][:,1,:], embedding_matrix[X_data][:,2,:]
                        , embedding_matrix[X_data][:,3,:], embedding_matrix[X_data][:,4,:], embedding_matrix[X_data][:,5,:]
                        , embedding_matrix[X_data][:,6,:], embedding_matrix[X_data][:,7,:], embedding_matrix[X_data][:,8,:]
                        , embedding_matrix[X_data][:,9,:], embedding_matrix[X_data][:,10,:], embedding_matrix[X_data][:,11,:]
                        , embedding_matrix[X_data][:,12,:], embedding_matrix[X_data][:,13,:], embedding_matrix[X_data][:,14,:]
                        , embedding_matrix[X_data][:,15,:], embedding_matrix[X_data][:,16,:], embedding_matrix[X_data][:,17,:]
                        , embedding_matrix[X_data][:,18,:], embedding_matrix[X_data][:,19,:], embedding_matrix[X_data][:,20,:]
                        , embedding_matrix[X_data][:,21,:], embedding_matrix[X_data][:,22,:], embedding_matrix[X_data][:,23,:]
                        , embedding_matrix[X_data][:,24,:], embedding_matrix[X_data][:,25,:], embedding_matrix[X_data][:,26,:]
                        , embedding_matrix[X_data][:,27,:]))  # shape (n, 840)

    return vector

