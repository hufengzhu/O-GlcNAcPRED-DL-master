import numpy as np
import math
import joblib
from sklearn import metrics
from keras.layers import *
from keras.models import *
from sklearn import preprocessing
from keras.callbacks import TensorBoard, EarlyStopping, ModelCheckpoint, ReduceLROnPlateau, Callback
from O_gly.AAindex.AAindex_sl import AAindex_hm_encoding, AAindex_ms_encoding, col_delete
from O_gly.word2vec.w2v_fea import w2v_fea_ms


scale = preprocessing.StandardScaler()

class LossHistory(Callback):
    def on_train_begin(self, logs={}):
        self.losses = []

    def on_batch_end(self, batch, logs={}):
        self.losses.append(logs.get('loss'))

def step_decay(epoch):
    initial_lrate = 0.0005
    drop = 0.8
    epochs_drop = 5.0
    lrate = initial_lrate * math.pow(drop, math.floor((1 + epoch) / epochs_drop))
    return lrate

# 1
def createModel1():

    word_input = Input(shape=(28, 29), name='word_input')

    overallResult = Convolution1D(filters=32, kernel_size=1, padding='same',activation="relu")(word_input)

    overallResult = Bidirectional(LSTM(50, dropout=0.5, activation='tanh', return_sequences=True))(overallResult)

    overallResult = Flatten()(overallResult)

    overallResult = Dense(32, activation='sigmoid')(overallResult)

    ss_output = Dense(1, activation='sigmoid', name='ss_output')(overallResult)

    return Model(inputs=[word_input], outputs=[ss_output])


# 3
def createModel3():

    word_input = Input(shape=(28, 30), name='word_input')

    overallResult = Convolution1D(filters=32, kernel_size=1, padding='same', activation="relu")(word_input)

    overallResult = Convolution1D(filters=16, kernel_size=1, padding='same', activation="relu", name='CNN_output')(overallResult)

    overallResult = Bidirectional(LSTM(50, dropout=0.5, activation='tanh', return_sequences=True),name='before_Dense')(overallResult)

    overallResult = Flatten()(overallResult)

    overallResult = Dense(32, activation='sigmoid')(overallResult)

    ss_output = Dense(1, activation='sigmoid', name='ss_output')(overallResult)

    return Model(inputs=[word_input], outputs=[ss_output])


# 4
def createModel4():

    word_input = Input(shape=(28, 29), name='word_input')

    overallResult = Convolution1D(filters=32, kernel_size=1, padding='same', activation="relu")(word_input)
    overallResult1 = Convolution1D(filters=32, kernel_size=1, padding='same', activation="relu")(word_input)

    overallResult = Bidirectional(LSTM(50, dropout=0.5, activation='tanh', return_sequences=True))(overallResult)
    overallResult1 = LSTM(32, return_sequences=True)(overallResult1)

    overallResult = Flatten()(overallResult)
    overallResult = Dropout(0.5)(overallResult)

    overallResult1 = Flatten()(overallResult1)
    overallResult1 = Dropout(0.5)(overallResult1)

    overallResult = Concatenate(axis=1)([overallResult, overallResult1])

    overallResult = Dense(64, activation='sigmoid')(overallResult)  # 之前是 128

    ss_output = Dense(1, activation='sigmoid', name='ss_output')(overallResult)

    return Model(inputs=[word_input], outputs=[ss_output])



seed = 7
batchSize = 32
maxEpochs = 200


## ============================================= 数据 ========================================
file = './test_data/independent test sets/Ind_M.fasta'


## 1-4 AAindex_gly
encodings = AAindex_ms_encoding(file)
X = np.array(encodings).reshape(-1, 841)
X = col_delete(X)
X1 = scale.fit_transform(X).reshape(X.shape[0], 28, 29)

## 5 w2v
w2v_mod = './feature/word2vec/w2v_ms_w4_v30.model'
X2 = w2v_fea_ms(file, w2v_mod).reshape(X.shape[0], 28, 30)


## ============================================= 模型 ========================================
M1 = createModel4()
M1.load_weights('./model/MS/model1.h5')
y_predict_test1 = M1.predict({'word_input': X1})

M2 = createModel1()
M2.load_weights('./model/MS/model2.h5')
y_predict_test2 = M2.predict({'word_input': X1})

M3 = createModel4()
M3.load_weights('./model/MS/model3.h5')
y_predict_test3 = M3.predict({'word_input': X1})

M4 = createModel4()
M4.load_weights('./model/MS/model4.h5')
y_predict_test4 = M4.predict({'word_input': X1})
#
M5 = createModel3()
M5.load_weights('./model/MS/model5.h5')
y_predict_test5 = M5.predict({'word_input': X2})


## ============================================= 集成 ========================================
sum_lst = []

for index, item in enumerate(y_predict_test1):

    sum_lst.append(((item*0.2) + (y_predict_test2[index])*0.05 + (y_predict_test3[index])*0.15 + (y_predict_test4[index])*0.3+ (y_predict_test5[index])*0.3))

b = np.array(sum_lst)

# 输出预测结果，>=0.5为O糖基化位点，<0.5不是
print(b)
