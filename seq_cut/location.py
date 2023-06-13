from Bio import SeqIO

for seq_record in SeqIO.parse("NUP62_cut.fasta", "fasta"):
    na = seq_record.id
    na_nam = na.split('|')[1]    # 肽段的名称
    na_pos = int(na.split('|')[-1].split('=')[-1])  # S/T的位置
    # se = seq_record.seq  # 序列

    print(na_nam)
    print(na_pos)
