# -*- coding: utf-8 -*-
"""
Pycharm Editor
Create by zhwei at 2020/03/18
Python：3.7.0
"""

import re
import argparse
import math
from collections import Counter
import pandas as pd
import numpy as np
from pathlib import Path
import sys

def read_Fasta(file: str):
    """
    :param file: fasta format
    :return: name and fasta sequences
    """
    #判断文件是否存在
    if Path(file).exists() == False:
        print(f'Error: {file} does not exist.')
        sys.exit(1)

    #读取文件内容，判断是否为fasta格式
    with open(file) as f:
        records = f.read()
    if re.search('>',records) == None:
        print(f"{file} seems not in fasta format")
        sys.exit(1)

    records = records.split('>')[1:]
    Fasta = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0],''.join(array[1:]).upper()

        #判断非自然氨基酸序列
        if re.search('[BJOUXZ]',sequence):
            print(f">{name} in {file} contains non-natural AA [BJOUXZ]")
            # sys.exit(1)

        #判断序列长度大于4
        elif len(sequence) <= 4:
            print(f"Sequence >{name} in {file} length <= 4")
            sys.exit(1)
        else:
            Fasta.append([name, sequence])
    return Fasta



def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range



def AAC(fastas):
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    header = []
    for i in AA:
        header.append(i)

    for sequence in fastas:
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = []
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return header, np.array(encodings, dtype=float)

def GAAC(fastas, **kw):
    group = {
        'alphatic': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharge': 'KRH',
        'negativecharge': 'DE',
        'uncharge': 'STCPNQ'
    }

    groupKey = group.keys()

    encodings = []
    header = []
    for key in groupKey:
        header.append(key)

    for sequence in fastas:
        code = []
        count = Counter(sequence)
        myDict = {}
        for key in groupKey:
            for aa in group[key]:
                myDict[key] = myDict.get(key, 0) + count[aa]

        for key in groupKey:
            code.append(myDict[key]/len(sequence))
        encodings.append(code)

    return header, np.array(encodings, dtype=float)

def GDPC(fastas, **kw):
    group = {
        'alphaticr': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharger': 'KRH',
        'negativecharger': 'DE',
        'uncharger': 'STCPNQ'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = [] + dipeptide

    for sequence in fastas:

        code = []
        myDict = {}
        for t in dipeptide:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 2 + 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                sequence[j + 1]]] + 1
            sum = sum + 1

        if sum == 0:
            for _ in dipeptide:
                code.append(0)
        else:
            for t in dipeptide:
                code.append(myDict[t] / sum)
        encodings.append(code)
    return header, np.array(encodings, dtype=float)

def Rvalue(aa1, aa2, AADict, Matrix):
    return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

def PAAC(fastas, lambdaValue=4, w=0.05, **kw):
    dataFile = './static/feature_data/PAAC.txt'

    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])

    encodings = []
    header = []
    for aa in AA:
        header.append('Xc1.' + aa)
    for n in range(1, lambdaValue + 1):
        header.append('Xc2.lambda' + str(n))

    for sequence in fastas:
        code = []
        theta = []
        for n in range(1, lambdaValue + 1):
            theta.append(
                sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (
                    len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
    return header, np.array(encodings, dtype=float)

def DPC(fastas, **kw):
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
    header = [] + diPeptides

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for sequence in fastas:
        code = []
        tmpCode = [0] * 400
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return header, np.array(encodings, dtype=float)

def CountCTDC(seq1, seq2):
    sum = 0
    for aa in seq1:
        sum = sum + seq2.count(aa)
    return sum

def CountCTDD(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence))
                    break
        if myCount == 0:
            code.append(0)
    return code



def CTD(fastas, **kw):
    group1 = {
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity':        'LIFWCMVY',
        'polarizability':  'GASDT',
        'charge':          'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess':   'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity':        'PATGS',
        'polarizability':  'CPNVEQIL',
        'charge':          'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess':   'RKQEND'
    }
    group3 = {
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity':        'HQRKNED',
        'polarizability':  'KMHFRYW',
        'charge':          'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess':   'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
    'hydrophobicity_FASG890101', 'normwaalsvolume',
    'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings_c = []
    encodings_t = []
    encodings_d = []
    headerc, headert, headerd = [], [], []
    for p in property:
        for g in range(1, len(groups) + 1):
            headerc.append(p + '.G' + str(g))
        for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
            headert.append(p + '.' + tr)
        for s in ('1', '2', '3'):
            for d in ['0', '25', '50', '75', '100']:
                headerd.append(p + '.' + s + '.residue' + d)

    for sequence in fastas:
        code_c = []
        code_t = []
        code_d = []
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        for p in property:
            c1 = CountCTDC(group1[p], sequence) / len(sequence)
            c2 = CountCTDC(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code_c = code_c + [c1, c2, c3]
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
            code_t = code_t + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]
            code_d = code_d + CountCTDD(group1[p], sequence) + CountCTDD(group2[p], sequence) + CountCTDD(group3[p], sequence)
        encodings_c.append(code_c)
        encodings_t.append(code_t)
        encodings_d.append(code_d)

    data_c = np.array(encodings_c)
    data_t = np.array(encodings_t)
    data_d = np.array(encodings_d)

    header = headerc + headert + headerd
    encodings = np.hstack((data_c, data_t, data_d))
    return header, encodings


def AAI(fastas):
    aaidx = pd.read_csv('./static/feature_data/aaindex1.txt', sep='\t')
    encodings = []
    header = list(aaidx.iloc[:,1])
    for i in fastas:
        s = 0
        code = []
        count = Counter(str(i))
        for jj in range(len(aaidx)):
            AA_idx = aaidx.iloc[jj,2:]
            for key in count:
                if key not in 'BJOUXZ':
                    s += count[key]*AA_idx[key]/len(str(i))
            code.append(s)
        encodings.append(code)
    return header, np.array(encodings, dtype=float)

def BPNC(fastas,name):
    letters = 'ACDEFGHIKLMNPQRSTVWY'

    aa_length = len(fastas[0])
    header = [f'BINARY.{name}.{i+1}' for i in range(aa_length*20)]

    encodings = []
    for fasta in fastas:
        encoding = []
        for AA in fasta:
            for AA1 in letters:
                tag = 1 if AA == AA1 else 0
                encoding.append(tag)
        encodings.append(encoding)
    return header, np.array(encodings, dtype=float)

def CT5(fastas):
    fastas_CT5  = [i[-5:] for i in fastas]
    return BPNC(fastas_CT5,'CT5')

def NT5(fastas):
    fastas_NT5  = [i[:5] for i in fastas]
    return BPNC(fastas_NT5,'NT5')


def all_feature(fastafile):
    temp = np.array(read_Fasta(fastafile))
    fastas_name, fastas = temp[:, 0], temp[:, 1]
    feature_list = ['AAI(fastas)','CT5(fastas)','CTD(fastas)','DPC(fastas)','GAAC(fastas)','GDPC(fastas)','NT5(fastas)','PAAC(fastas)']

    #对header 和 encodings初始化
    header, encodings = AAC(fastas)
    subfeature_index = [encodings.shape[-1]]
    for i in feature_list:
        sub_header, sub_feature = eval(i)
        header += sub_header
        encodings = np.hstack((encodings, sub_feature))
        subfeature_index.append(sub_feature.shape[-1])

    # 子特征索引
    feature_index = []
    idx = 0
    for i in subfeature_index:
        feature_index.append([m + idx for m in range(i)])
        idx += i

    return pd.DataFrame(encodings,index=fastas_name,columns=header), feature_index


def get_args():
    """
    Usage:
    python Features.py -f file.fasta -o ./data/Features/Base.csv
    """

    parser = argparse.ArgumentParser(usage="Usage Tip;",
                                     description="APPred Feature Extraction")
    parser.add_argument("--file", "-f", help="input training file(.fasta)")
    parser.add_argument("--out", "-o", help="output path and filename")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    features,subfeature_index = all_feature(args.file)
    np.savetxt('./data/Features/fea_idx.txt',subfeature_index,delimiter=',',fmt='%d')
    features.to_csv(args.out)
    print(subfeature_index)