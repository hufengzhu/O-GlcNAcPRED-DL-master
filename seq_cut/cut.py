from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# 截取长为29的肽段
def strcat(seq, front=True, default_char='X', default_len=29):
    temp = default_char * (default_len - len(seq))
    if front:
        return temp + seq
    else:
        return seq + temp


def na_se(file_name):
    result = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        na = seq_record.id
        se = seq_record.seq
        if len(se) < 29:
            print(na)
        for i in range(len(str(se))):
            if se[i] == 'S' or se[i] == 'T':
                if i>= 14 and i+14 < len(se):
                    new_record = SeqRecord(se[i-14:i+15], id=na+f"|position={i}")
                elif i < 14:
                    new_record = SeqRecord(strcat(se[:i+15]), id=na + f"|position={i}")
                elif i+14 >= len(se):
                    new_record = SeqRecord(strcat(se[i-14:], False), id=na + f"|position={i}")
                result.append(new_record)
    with open("NUP62_cut.fasta", "w") as fasta_file:   # 截取后的数据
        SeqIO.write(result, fasta_file, "fasta")

    return fasta_file

na_se('NUP62.fasta')  # 要截取的数据

