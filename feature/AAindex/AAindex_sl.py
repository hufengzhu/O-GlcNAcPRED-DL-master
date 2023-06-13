import numpy as np
import re
from Bio import SeqIO


def AAindex_hm_encoding(file_name):
    with open(r'D:\python code\O-GlcNAcPRED-DL-master\feature\AAindex\AAindex_normalized.txt') as f:
        records = f.readlines()[1:]
    AA_aaindex = 'ARNDCQEGHILKMFPSTWYVX'
    AAindex = []
    AAindexName = []
    for i in records:
        AAindex.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
        AAindexName.append(i.rstrip().split()[0] if i.rstrip() != '' else None)
    props = 'HUTJ700101:VHEG790101:ISOY800106:ROBB760113:CIDH920101:OOBM770102:TANS770102:ROSM880101:QIAN880113:SNEP660104:OOBM850103:MAXF760106:BAEK050101:RADA880101:HOPA770101:FAUJ880112:GEIM800101:WOLR810101:CHOP780207:FASG760105:FAUJ830101:CHAM810101:QIAN880124:BASU050102:GRAR740101:AURR980116:AURR980117:WOLS870101:ZIMJ680104'.split(
        ':')  ## human挑选出的29个理化特征
    if props:
        tmpIndexNames = []
        tmpIndex = []
        for p in props:
            if AAindexName.index(p) != -1:
                tmpIndexNames.append(p)
                tmpIndex.append(AAindex[AAindexName.index(p)])
        if len(tmpIndexNames) != 0:
            AAindexName = tmpIndexNames
            AAindex = tmpIndex

    index = {}
    for i in range(len(AA_aaindex)):
        index[AA_aaindex[i]] = i

    encoding_aaindex = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        sequence = seq_record.seq
        code = []
        for aa in sequence:
            if aa == 'X' or aa == 'U':
                for j in AAindex:
                    code.append(0)
                continue
            for j in AAindex:
                code.append(j[index[aa]])
        encoding_aaindex.append(code)
    return encoding_aaindex


def AAindex_ms_encoding(file_name):
    with open(r'D:\python code\O-GlcNAcPRED-DL-master\feature\AAindex\AAindex_normalized.txt') as f:
        records = f.readlines()[1:]
    AA_aaindex = 'ARNDCQEGHILKMFPSTWYVX'
    AAindex = []
    AAindexName = []
    for i in records:
        AAindex.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
        AAindexName.append(i.rstrip().split()[0] if i.rstrip() != '' else None)
    props = 'BAEK050101:FAUJ830101:LEVM760105:LAWE840101:QIAN880134:FAUJ880104:WERD780101:GEIM800102:TANS770102:BEGF750103:KRIW790101:HOPA770101:CHAM830104:PONP800108:EISD860102:BASU050102:LEVM760103:QIAN880124:GEOR030107:ZIMJ680104:AURR980117:CHAM810101:MITS020101:KARP850102:CHOP780207:HUTJ700101:LEVM760106:JOND750102:RICJ880110'.split(
        ':')  ## mouse挑选出的29个理化特征
    if props:
        tmpIndexNames = []
        tmpIndex = []
        for p in props:
            if AAindexName.index(p) != -1:
                tmpIndexNames.append(p)
                tmpIndex.append(AAindex[AAindexName.index(p)])
        if len(tmpIndexNames) != 0:
            AAindexName = tmpIndexNames
            AAindex = tmpIndex

    index = {}
    for i in range(len(AA_aaindex)):
        index[AA_aaindex[i]] = i

    encoding_aaindex = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        sequence = seq_record.seq
        code = []
        for aa in sequence:
            if aa == 'X':
                for j in AAindex:
                    code.append(0)
                continue
            for j in AAindex:
                code.append(j[index[aa]])
        encoding_aaindex.append(code)
    return encoding_aaindex


def col_delete(data):  # delete the columns woth same elements
    col_del = [406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426,
               427, 428, 429, 430, 431, 432, 433, 434]
    data_new = np.delete(data, col_del, axis=1)
    return data_new

