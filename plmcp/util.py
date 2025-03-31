import os
from typing import Union, List, Tuple, Dict

import torch
import numba
from tqdm import tqdm

import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import re
#huyue
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_fasta(fn_fasta):
    prot2seq = {}
    with open(fn_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)
            prot = record.id
            prot2seq[prot] = seq
    return list(prot2seq.keys()), prot2seq
def read_CP_fasta(fn_fasta):
    prot2seq = {}
    listseq = []
    listprot = []
    listrecord = []
    with open(fn_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)
            seq= seq + seq #huyue 
            prot = record.id
            prot2seq[prot] = seq
            listseq.append(seq)
            listprot.append(prot)
            record.seq=record.seq+record.seq
            listrecord.append(record)
        #short_seq_iterator = (record for record in SeqIO.parse(handle, "fasta"): record.seq=record.seq+record.seq)
    #input_seq_iterator = SeqIO.parse("", "genbank")
    #short_seq_iterator = (record for record in SeqIO.parse(handle, "fasta") record.seq=record.seq+record.seq)
    #SeqIO.write(short_seq_iterator, "short_seqs.fasta", "fasta")
    #record = SeqRecord(Seq(seq), id=prot, description="CPsequence")
    #record = SeqRecord(Seq(seq), id=prot, description="")
    #record = SeqRecord(Seq(listseq), id=listprot, description="")
    #output_file = "example/CP_output.fasta"
    #SeqIO.write(record, output_file, "fasta")
    SeqIO.write(listrecord, "example/CP_output.fasta", "fasta")

    return list(prot2seq.keys()), prot2seq
def read_CP_target_fasta(fn_fasta):
    prot2seq = {}
    listrecord = []
    with open(fn_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)
            seq= seq + seq #huyue 
            prot = record.id
            prot2seq[prot] = seq
            record.seq=record.seq+record.seq
            listrecord.append(record)
    #record = SeqRecord(Seq(seq), id=prot, description="CPsequence")
    #record = SeqRecord(Seq(seq), id=prot, description="")
    #record = SeqRecord(Seq(prot2seq), id=prot, description="")
    #output_file = "example/CP_target_output.fasta"
    #SeqIO.write(record, output_file, "fasta")
    SeqIO.write(listrecord, "example/CP_target_output.fasta", "fasta")

    return list(prot2seq.keys()), prot2seq

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO

def read_target_fasta_to_files(fn_fasta):
    prot2seq = {}
    listrecord_target = []
    listrecord_CP = []
    # 如果输入是一行字符串，将其转换为伪装的 FASTA 格式
    if not fn_fasta.startswith(">"):
        fn_fasta = ">sequenceTarget\n" + fn_fasta       
    # Use StringIO to create a file-like object from the input string
    handle = StringIO(fn_fasta)
    # Parse the input FASTA string
    seq_list = []  # Default value to avoid UnboundLocalError
    seq_CP_list = []  # Default value to avoid UnboundLocalError    
    # Parse the input FASTA string

	


    # 解析输入的 FASTA 字符串
    for record in SeqIO.parse(handle, "fasta"):
        # 提取原始序列并去除首尾空白字符
        seq = str(record.seq).strip()
        seq = re.sub(r'\s+', '', seq)
        seq_list.append(seq)

        # 保存去除空白字符后的原始序列
        prot = record.id
        prot2seq[prot] = seq
        record_target = SeqRecord(Seq(seq), id=prot, description=record.description)
        listrecord_target.append(record_target)

        # 创建一个新的序列，将原始序列复制两次
        seq_CP = seq + seq  # 复制序列
        seq_CP_list.append(seq_CP)
        
        # 创建新的 SeqRecord 对象，包含复制后的序列
        record_CP = SeqRecord(Seq(seq_CP), id=prot, description="TargetCPsequence")
        listrecord_CP.append(record_CP)	
    # Write the sequences to two separate FASTA files
    SeqIO.write(listrecord_target, "example/target_output.fasta", "fasta")
    SeqIO.write(listrecord_CP, "example/CP_target_output.fasta", "fasta")	
	
	
    # Return the list of protein IDs and the dictionary mapping IDs to sequences
    return seq_list, seq_CP_list

def read_query_fasta_to_files(fn_fasta):
    prot2seq = {}
    listrecord_target = []
    listrecord_CP = []
    # 如果输入是一行字符串，将其转换为伪装的 FASTA 格式
    if not fn_fasta.startswith(">"):
        fn_fasta = ">sequenceQuery\n" + fn_fasta    
    # Use StringIO to create a file-like object from the input string
    handle = StringIO(fn_fasta)
    # Parse the input FASTA string
    seq_list = []  # Default value to avoid UnboundLocalError
    seq_CP_list = []  # Default value to avoid UnboundLocalError      

	
    # 解析输入的 FASTA 字符串
    for record in SeqIO.parse(handle, "fasta"):
        # 提取原始序列并去除首尾空白字符
        seq = str(record.seq).strip()
        seq = re.sub(r'\s+', '', seq)
        seq_list.append(seq)

        # 保存去除空白字符后的原始序列
        prot = record.id
        prot2seq[prot] = seq
        record_target = SeqRecord(Seq(seq), id=prot, description=record.description)
        listrecord_target.append(record_target)

        # 创建一个新的序列，将原始序列复制两次
        seq_CP = seq + seq  # 复制序列
        seq_CP_list.append(seq_CP)
        
        # 创建新的 SeqRecord 对象，包含复制后的序列
        record_CP = SeqRecord(Seq(seq_CP), id=prot, description="QueryCPsequence")
        listrecord_CP.append(record_CP)			
		
		
    # Write the sequences to two separate FASTA files
    SeqIO.write(listrecord_target, "example/query_output.fasta", "fasta")
    SeqIO.write(listrecord_CP, "example/CP_query_output.fasta", "fasta")

    # Return the list of protein IDs and the dictionary mapping IDs to sequences
    return seq_list, seq_CP_list
	
def make_parent_dir(path):
    filepath = Path(path)
    filepath.parent.mkdir(parents=True, exist_ok=True)

@numba.njit('f4[:,:](f4[:,:], f4[:,:])', nogil=True, fastmath=True, cache=True)
def dot_product(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
	assert X.ndim == 2 and Y.ndim == 2
	assert X.shape[1] == Y.shape[1]

	xlen: int = X.shape[0]
	ylen: int = Y.shape[0]
	embdim: int = X.shape[1]

	emb1_normed: np.ndarray = np.ones((xlen, embdim), dtype=np.float32)
	emb2_normed: np.ndarray = np.ones((ylen, embdim), dtype=np.float32)
	density: np.ndarray = np.empty((xlen, ylen), dtype=np.float32)
	# numba does not support sum() args other then first
	emb1_normed = X / 1
	emb2_normed = Y / 1
	density = emb1_normed @ emb2_normed.T
	return density

@numba.njit('f4[:,:](f4[:,:], f4[:,:])', nogil=True, fastmath=True, cache=True)
def embedding_cos_similarity(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
	assert X.ndim == 2 and Y.ndim == 2
	assert X.shape[1] == Y.shape[1]

	xlen: int = X.shape[0]
	ylen: int = Y.shape[0]
	embdim: int = X.shape[1]
	# normalize
	emb1_norm: np.ndarray = np.empty((xlen, 1), dtype=np.float32)
	emb2_norm: np.ndarray = np.empty((ylen, 1), dtype=np.float32)
	emb1_normed: np.ndarray = np.empty((xlen, embdim), dtype=np.float32)
	emb2_normed: np.ndarray = np.empty((ylen, embdim), dtype=np.float32)
	density: np.ndarray = np.empty((xlen, ylen), dtype=np.float32)
	# numba does not support sum() args other then first
	emb1_norm = np.expand_dims(np.sqrt(np.power(X, 2).sum(1)), 1)
	emb2_norm = np.expand_dims(np.sqrt(np.power(Y, 2).sum(1)), 1)
	emb1_normed = X / emb1_norm
	emb2_normed = Y / emb2_norm
	density = emb1_normed @ emb2_normed.T
	return density

def get_prefilter_list(prefilter_result, query_num):
    prefilter_list = []
    got_query = set()
    with open(prefilter_result) as fp:
        for line in tqdm(fp, desc='Load Result'):
            line_list = line.strip().split('\t')
            protein1 = line_list[0].split('.pdb')[0]
            got_query.add(protein1)
            if (len(got_query) > query_num):
                break
            protein2 = line_list[1].split('.pdb')[0]
            score = eval(line_list[2])
            prefilter_list.append(((protein1, protein2), score))
    return prefilter_list

def filter_result_dataframe(data: pd.DataFrame,
							column: Union[str, List[str]] = ['score']) -> \
								pd.DataFrame:
	'''
	keep spans with biggest score and len
	Args:
		data: (pd.DataFrame)
	Returns:
		filtred frame sorted by score
	'''
	data = data.sort_values(by=['len'], ascending=False)
	indices = data.indices.tolist()
	data['y1'] = [yx[0][0] for yx in indices]
	data['x1'] = [yx[0][1] for yx in indices]
	data['score'] = data['score'].round(3)

	if isinstance(column, str):
		column = [column]
	resultsflt = list()
	iterator = data.groupby(['y1', 'x1'])
	for col in column:
		for groupid, group in iterator:
			tmp = group.nlargest(1, [col], keep='first')
			resultsflt.append(tmp)
	resultsflt = pd.concat(resultsflt)
	# drop duplicates sometimes
	resultsflt = resultsflt.drop_duplicates(
		subset=['pathid', 'i', 'len', 'score'])
	# filter
	resultsflt = resultsflt.sort_values(by=['score'], ascending=False)
	return resultsflt

def embedding_load(fasta, embedding_path):
	_, sequences = read_fasta(fasta)
	embedding_result_dic = {}
	for single_sequence in sequences:
		embedding_path = Path(embedding_path)
		embedding_file = embedding_path / f"{single_sequence}.pt"
		embedding_result_dic[single_sequence] = torch.load(embedding_file)
	return embedding_result_dic

def draw_alignment(coords: List[Tuple[int, int]], seq1: str, seq2: str, output: Union[None, str]) -> str:
	'''
	draws alignment based on input coordinates
	Args:
		coords: (list) result of align list of tuple indices
		seq1: (str) full residue sequence 
		seq2: (str) full residue sequence
		output: (str or bool) if None output is printed
	'''
	assert isinstance(seq1, str) or isinstance(seq1[0], str), 'seq1 must be sting like type'
	assert isinstance(seq2, str)or isinstance(seq1[0], str), 'seq2 must be string like type'
	assert len(seq1) > 1 and len(seq2), 'seq1 or seq1 is too short'

	# check whether alignment indices exeed sequence len
	last_position = coords[-1]
	lp1, lp2 = last_position[0], last_position[1]
	if lp1 >= len(seq1):
		raise KeyError(f'mismatch between seq1 length and coords {lp1} - {len(seq1)} for seq2 {lp2} - {len(seq2)}')
	if lp2 >= len(seq2):
		raise KeyError(f'mismatch between seq1 length and coords {lp2} - {len(seq2)}')

	# container
	#alignment = dict(up=[], relation=[], down=[],upIndex=[],downIndex=[])
	#alignment = dict(up=[], relation=[], down=[], Index=[])
	alignment = dict(up=[], relation=[], down=[])
	c1_prev, c2_prev = -1, -1
	
	for c1, c2 in coords:
		# check if gap occur
		up_increment   = True if c1 != c1_prev else False
		down_increment = True if c2 != c2_prev else False
		
		if up_increment:
			up = seq1[c1]
		else:
			up = '-'

		if down_increment:
			down = seq2[c2]
		else:
			down = '-'

		if up_increment and down_increment:
			relation = '|'
		else:
			relation = ' '
			
		alignment['up'].append(up)
		alignment['relation'].append(relation)
		alignment['down'].append(down)
                #indexC1C2 = str(c1)
		#indexC1C2 = c1 + '\t'
                #indexC1C2=str(c1) + '\t' + str(c2) + ';'
		#alignment['Index'].append(indexC1C2)
                #alignment['Index'].append(indexC1C2)
		#string = ''.join(alignment['up']) + '\n'
                #alignment['Index'].append(c1','c2';')
                #alignment['Index'].append(c1)
		#alignment['downIndex'].append(c2)
			
		c1_prev = c1
		c2_prev = c2
	# merge into 3 line string
	if output != 'html':
		string = ''.join(alignment['up']) + '\n'
		string += ''.join(alignment['relation']) + '\n'
		string += ''.join(alignment['down'])
		#string += ''.join(alignment['Index'])
		if output is not None:
			return string
		else:
			print(string)
