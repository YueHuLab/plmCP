import itertools

import pandas as pd
from typing import Union
import numpy as np
from tqdm import tqdm

from .plmcp_util.alignment import plmcp_gather_all_paths, plmcp_search_paths
from .util import filter_result_dataframe, dot_product, read_fasta, embedding_load, embedding_cos_similarity, make_parent_dir, draw_alignment, read_CP_fasta,read_CP_target_fasta
from .embedding_generate import esm_embedding_generate

import scipy
class plmcp:
    '''
    main class for handling alignment extaction
    '''
    NORM: Union[bool, str] = True
    #NORM: Union[bool, str] = False
    MODE: str = 'local'
    #GAP_EXT: float = 1.0
    GAP_EXT: np.float32 = 1.0 #huyue

    def __init__(self, *args, **kw_args):
        pass

    def embedding_to_span(self, X, Y, result_mode : str = 'results') -> pd.DataFrame:

        ### PLMAlign (dot)
        X = X.numpy()
        Y = Y.numpy()
        densitymap = dot_product(X, Y)

        densitymap = densitymap.T
        
        #path = plmcp_gather_all_paths(densitymap,
        path,scorematrix = plmcp_gather_all_paths(densitymap,
                                norm=self.NORM,
                                mode=self.MODE,
                                gap_extension=self.GAP_EXT,
                                #with_scores = True if result_mode == 'all' else False)
                                with_scores = True)
        if result_mode == 'all':
            scorematrix = path[1]
            path = path[0]
            #huyue
        #results = plmcp_search_paths(densitymap
        #print(f"path \t path = {path}\n")
        results = plmcp_search_paths(scorematrix,
                               path=path,
                               mode=self.MODE,
                               #gap_extension=self.GAP_EXT,
                               as_df=True)
        if result_mode == 'all':
            return (results, densitymap, path, scorematrix)
        else:
            return results


def pairwise_align(embedding1, embedding2, seq1, seq2, mode, method = 'plmcp'):
    if method == 'plmcp':
        extr = plmcp()
        extr.MODE = mode
       # result_model='all'
        #results = extr.embedding_to_span(embedding2, embedding1, result_model)
        results = extr.embedding_to_span(embedding2, embedding1)
    else:
        assert method in {"plmcp"}
    #huyue
    df=results
    #df=df.sort_values(by='score_start',ascending=False)
    #df=df.sort_values(by='score',ascending=False)
    # Print best alignment
    #row = results.iloc[0]
    row = df.iloc[0]
    #row = df.iloc[1] #huyue

    aln = draw_alignment(row.indices, seq1, seq2, output='str')

    #print("row indeices\n")
    #print(row.indices)
    #print("\n")
    #print("seq1\n")
    #print(seq1)
    #print("\n")
    #print("seq2\n")
    #print(seq2)
    #print("\n")
    return row['score'].item(), aln,  row['alnlen'].item(),row.indices,row['scoreAll']

def plmCP_fisher_pipeline(query_fasta, target_fasta, mode = 'local', query_embedding_path = None, target_embedding_path = None, search_result_setting = None, output_path = None, if_stdout = True):
#def plmCP_fisher_pipeline(query_fasta, target_fasta, mode = 'global', query_embedding_path = None, target_embedding_path = None, search_result_setting = None, output_path = None, if_stdout = True):
    method = 'plmcp'
    print(f"Align with method: {method}")

    _, query_sequences = read_fasta(query_fasta)
    _, target_sequences = read_fasta(target_fasta)
    _, query_CP_sequences = read_CP_fasta(query_fasta)
    _, target_CP_sequences = read_CP_target_fasta(target_fasta)
    #
    #length_query=len(query_sequences[0])
    #print(f"length query{length_query}\n")
    #
    query_CP_fasta= "example/CP_output.fasta"
    target_CP_fasta= "example/CP_target_output.fasta"
    if (query_embedding_path == None):
        query_embeddings = esm_embedding_generate(query_fasta)
    else:
        query_embeddings = embedding_load(query_fasta, query_embedding_path)
    #query CP huyue
    query_CP_embeddings = esm_embedding_generate(query_CP_fasta)
    target_CP_embeddings = esm_embedding_generate(target_CP_fasta)

    if (target_embedding_path == None):
        target_embeddings = esm_embedding_generate(target_fasta)
    else:
        target_embeddings = embedding_load(target_fasta, target_embedding_path)
    #print 
    #length_query=list(query_embeddings.values())[0].shape()
    length_query= []
    length_target= []
    for i in range(len(list(query_embeddings.values()))):
            length_query.append(len(list(query_embeddings.values())[i]))
    for i in range(len(list(target_embeddings.values()))):
            length_target.append(len(list(target_embeddings.values())[i]))
    #length_query=len(list(query_embeddings.values())[0])
    #length_target=len(list(target_embeddings.values())[0])
    #print(f"length query{length_query}\n")
    #print(f"length target{length_target}\n")
    ll_query=len(length_query)
    ll_target=len(length_target)
    new_length_query=[i for i in length_query for j in range(ll_target)]
    new_length_target=length_target*ll_query



    if (output_path != None):
        output_score = output_path + 'score'
        make_parent_dir(output_score)
        output_alignment = output_path + 'alignment'
        f1 = open(output_score, 'w')
        f2 = open(output_alignment, 'w')
        
        protein_pair_dict = {}
        for protein in query_sequences:
            protein_pair_dict[protein] = []
        output_score_sort = output_path + 'score_sort'
    LISTscore= []
    LISTresults= []
    LISTalnlen= []
    LISTrowIndex= []
    LISTscoreAll= []
    LISTscore_tCP1= []
    LISTresults_tCP1= []
    LISTalnlen_tCP1= []
    LISTrowIndex_tCP1= []
    LISTscoreAll_tCP1= []

    LISTquery_seq=[]
    LISTtarget_seq=[]
    LISTquery=[]
    LISTtarget=[]

    if (search_result_setting == None):
        for single_query in tqdm(query_sequences, desc="Query"):
            for single_target in target_sequences:
                score, results,alnlen,rowIndex,scoreAll = pairwise_align(query_embeddings[single_query], target_embeddings[single_target], query_sequences[single_query], target_sequences[single_target], mode, method = method)
                LISTscore.append(score)
                LISTresults.append(results)
                LISTalnlen.append(alnlen)
                LISTrowIndex.append(rowIndex)
                LISTscoreAll.append(scoreAll)
                LISTquery.append(single_query)
                LISTtarget.append(single_target)
                LISTquery_seq.append(query_sequences[single_query])
                LISTtarget_seq.append(target_sequences[single_target])
                #print(f"query \t query {single_query}\n")
                if if_stdout:
                    print(f"{single_query}\t{single_target}\t Score = {score}\n")
                    print(f"{single_query}\t{single_target}\n{results}\n")
                    print(f"{single_query}\t{single_target}\t Len={alnlen}\n")
                    #print(f"{single_query}\t{single_target}\t ScoreAll = {scoreAll}\n")
                if (output_path != None):
                    f1.write(f"{single_query}\t{single_target}\t{score}\n")
                    f2.write(f"{single_query}\t{single_target}\n{results}\n\n")
                    protein_pair_dict[single_query].append((single_target, score))
            for single_CP_target in target_CP_sequences:
                score_tCP1, results_tCP1,alnlen_tCP1,rowIndex_tCP1,scoreAll_tCP1 = pairwise_align(query_embeddings[single_query], target_CP_embeddings[single_CP_target], query_sequences[single_query], target_CP_sequences[single_CP_target], mode, method = method)
                LISTscore_tCP1.append(score_tCP1)
                LISTresults_tCP1.append(results_tCP1)
                LISTalnlen_tCP1.append(alnlen_tCP1)
                LISTrowIndex_tCP1.append(rowIndex_tCP1)
                LISTscoreAll_tCP1.append(scoreAll_tCP1)

    else:
        search_result_path = search_result_setting[0]
        top = search_result_setting[1]

        with open(search_result_path, "r") as f:
            pairs = f.readlines()

        for line in tqdm(pairs, desc="Search result"):
            single_query, single_target, similarity = line.strip().split()
            similarity = eval(similarity)

            if ((top != None) and (similarity<top)):
                continue

            score, results = pairwise_align(query_embeddings[single_query], target_embeddings[single_target], query_sequences[single_query], target_sequences[single_target], mode, method = method)
            if if_stdout:
                print(f"{single_query}\t{single_target}\t Score = {score}\n")
                print(f"{single_query}\t{single_target}\n{results}\n")
            if (output_path != None):
                f1.write(f"{single_query}\t{single_target}\t{score}\n")
                f2.write(f"{single_query}\t{single_target}\n{results}\n\n")
                protein_pair_dict[single_query].append((single_target, score))
    
    if (output_path != None):
        f1.close()
        f2.close()
    
        for query_protein in query_sequences:
            protein_pair_dict[query_protein] = sorted(protein_pair_dict[query_protein], key=lambda x:x[1], reverse=True)
        
        with open(output_score_sort, 'w') as f3:
            for query_protein in query_sequences:
                for pair in protein_pair_dict[query_protein]:
                    f3.write(f"{query_protein}\t{pair[0]}\t{pair[1]}\n")

########################################
    if (output_path != None):
        output_CP_score = output_path + 'CPscore'
        make_parent_dir(output_CP_score)
        output_CP_alignment = output_path + 'CPalignment'
        fCP1 = open(output_CP_score, 'w')
        fCP2 = open(output_CP_alignment, 'w')
        
        protein_CP_pair_dict = {}
        for CPprotein in query_CP_sequences:
            protein_CP_pair_dict[CPprotein] = []
        output_CP_score_sort = output_path + 'CPscore_sort'
    LISTCPscore= []
    LISTCPresults= []
    LISTCPalnlen= []
    LISTCProwIndex= []
    LISTCPscoreAll= []
    LISTCPscore_tCP2= []
    LISTCPresults_tCP2= []
    LISTCPalnlen_tCP2= []
    LISTCProwIndex_tCP2= []
    LISTCPscoreAll_tCP2= []

    if (search_result_setting == None):
        for single_CP_query in tqdm(query_CP_sequences, desc="Query"):
            for single_target in target_sequences:
                CPscore, CPresults,CPalnlen,rowIndexCP,scoreAllCP = pairwise_align(query_CP_embeddings[single_CP_query], target_embeddings[single_target], query_CP_sequences[single_CP_query], target_sequences[single_target], mode, method = method)
                LISTCPscore.append(CPscore)
                LISTCPresults.append(CPresults)
                LISTCPalnlen.append(CPalnlen)
                LISTCProwIndex.append(rowIndexCP)
                LISTCPscoreAll.append(scoreAllCP)
                if if_stdout:
                    print(f"{single_CP_query}\t{single_target}\t CpScore = {CPscore}\n")
                    print(f"{single_CP_query}\t{single_target}\n cpResult = {CPresults}\n")
                    print(f"{single_CP_query}\t{single_target}\t cpLen = {CPalnlen}\n")
                    #print(f"{single_CP_query}\t{single_target}\t ScoreAll ={scoreAllCP}\n")
                if (output_path != None):
                    fCP1.write(f"{single_CP_query}\t{single_target}\t{CPscore}\n")
                    fCP2.write(f"{single_CP_query}\t{single_target}\n{CPresults}\n\n")
                    protein_CP_pair_dict[single_CP_query].append((single_target, CPscore))
            for single_CP_target in target_CP_sequences:
                CPscore_tCP2, CPresults_tCP2,CPalnlen_tCP2,CProwIndex_tCP2,CPscoreAll_tCP2 = pairwise_align(query_CP_embeddings[single_CP_query], target_CP_embeddings[single_CP_target], query_CP_sequences[single_CP_query], target_CP_sequences[single_CP_target], mode, method = method)
                LISTCPscore_tCP2.append(CPscore_tCP2)
                LISTCPresults_tCP2.append(CPresults_tCP2)
                LISTCPalnlen_tCP2.append(CPalnlen_tCP2)
                LISTCProwIndex_tCP2.append(CProwIndex_tCP2)
                LISTCPscoreAll_tCP2.append(CPscoreAll_tCP2)

    else:
        search_result_path = search_result_setting[0]
        top = search_result_setting[1]
#huyue
        with open(search_result_path, "r") as fCP:
            CPpairs = fCP.readlines()

        for CPline in tqdm(CPpairs, desc="Search result"):
            single_CP_query, single_target, CPsimilarity = CPline.strip().split()
            CPsimilarity = eval(CPsimilarity)

            if ((top != None) and (CPsimilarity<top)):
                continue

            CPscore, CPresults = pairwise_align(query_CP_embeddings[single_CP_query], target_embeddings[single_target], query_CP_sequences[single_CP_query], target_sequences[single_target], mode, method = method)
            if if_stdout:
                print(f"{single_CP_query}\t{single_target}\t CPScore = {CPscore}\n")
                print(f"{single_CP_query}\t{single_target}\n{CPresults}\n")
            if (output_path != None):
                fCP1.write(f"{single_CP_query}\t{single_target}\t{CPscore}\n")
                fCP2.write(f"{single_CP_query}\t{single_target}\n{CPresults}\n\n")
                protein_CP_pair_dict[single_CP_query].append((single_target, CPscore))
    
    if (output_path != None):
        fCP1.close()
        fCP2.close()
    
        for query_CP_protein in query_CP_sequences:
            protein_CP_pair_dict[query_CP_protein] = sorted(protein_CP_pair_dict[query_CP_protein], key=lambda x:x[1], reverse=True)
        
        with open(output_CP_score_sort, 'w') as fCP3:
            for query_CP_protein in query_CP_sequences:
                for CPpair in protein_CP_pair_dict[query_CP_protein]:
                    fCP3.write(f"{query_CP_protein}\t{CPpair[0]}\t{CPpair[1]}\n")

########################################
    #if (score<CPscore):
    #     print(f"The two proteins maybe evolved by circular permutation\n")

    #total_length=length_query+length_target
    #n1=total_length-alnlen
    #A=alnlen_tCP1- 2*length_target
    #B=CPalnlen-2*length_query
    #CPcheck=(A*B)/(length_query*length_target)

    assert len(LISTrowIndex_tCP1)== len(LISTCProwIndex), 'two list are not equl'
    #assert len(LISTrowIndex_tCP1)== len(new_length_query), 'two list are not equl'
    #assert len(new_length_query)== len(new_length_target), 'two list are not equl'
    for i in range(len(LISTrowIndex_tCP1)):
    #for indexTcp1 in LISTrowIndex_tCP1:
    #    for indexCP in LISTCProwIndex:
        c1_prev, c2_prev = -1, -1
        Arelation=0
        indexTcp1=LISTrowIndex_tCP1[i]
        indexCP=LISTCProwIndex[i]
        for c1, c2 in indexTcp1:
        # check if gap occur
            up_increment   = True if c1 != c1_prev else False
            down_increment   = True if c1 != c1_prev else False
            if up_increment and down_increment:
                Arelation =Arelation + 1
            c1_prev = c1
            c2_prev = c2
	# merge into 3 line string

        c1_prev, c2_prev = -1, -1
        Brelation=0	
        for c1, c2 in indexCP:
        # check if gap occur
            up_increment   = True if c1 != c1_prev else False
            down_increment   = True if c1 != c1_prev else False
            if up_increment and down_increment:
                Brelation = Brelation+1
            c1_prev = c1
            c2_prev = c2
	# merge into 3 line string
    

    #A=alnlen_tCP1-Arelation
    #B=CPalnlen-Brelation
#        print(f"Test between  {LISTquery[i]} \t {LISTtarget[i]}\n")
        A=Arelation
        B=Brelation
        #CPcheck=(A*B)/(length_query*length_target)
        CPcheck=(A*B)/(len(LISTquery_seq[i])*len(LISTtarget_seq[i]))


 #       print(f"resFIsher CPcheck {CPcheck}\n")
        #print(f"length query and target {new_length_query[i]}\t{new_length_target[i]}\n")
 #       print(f"length query and target {len(LISTquery_seq[i])}\t{len(LISTtarget_seq[i])}\n")
 #       if (score<CPscore and alnlen < alnlen_tCP1 and CPcheck > 0.5):
 #           print(f"The two proteins maybe evolved by circular permutation\n")


        #CPcheck2=(1-A/length_query)*(1-B/length_target)
        #CPcheck2=(1-A/new_length_query[i])*(1-B/new_length_target[i])
        CPcheck2=(1-A/len(LISTquery_seq[i]))*(1-B/len(LISTtarget_seq[i]))
        #LISTscore[i]
        #LISTresults.append(results)
        #LISTalnlen[i]
        
        #LISTscore_tCP1[i]
        #LISTresults_tCP1.append(results_tCP1)
        #LISTalnlen_tCP1[i]
        
        #LISTCPscore.append[i]
        #LISTCPresults.append(CPresults)
        #LISTCPalnlen.append[i]

        #LISTCPscore_tCP2.append[i]
        #LISTCPresults_tCP2.append(CPresults_tCP2)
        #LISTCPalnlen_tCP2.append[i]


#        print(f"resFIsher CPcheck2 {CPcheck2}\n")
        print(f"TestBetween \t {LISTquery[i]} \t {LISTtarget[i]}\t queryLen \t{len(LISTquery_seq[i])}\t targetLen \t{len(LISTtarget_seq[i])}\t CP1 \t {CPcheck}\t CP2 \t{CPcheck2}\t info \t { LISTscore[i]}\t {LISTalnlen[i]}\t {LISTscore_tCP1[i]} \t {LISTalnlen_tCP1[i]} \t {LISTCPscore[i]} \t {LISTCPalnlen[i]} \t {LISTCPscore_tCP2[i]} \t {LISTCPalnlen_tCP2[i]} \n")
#        if (score<CPscore and alnlen < alnlen_tCP1 and CPcheck2 < 0.5):
#            print(f"The two proteins maybe evolved by circular permutation\n,check2")
def plmCP_mutil_pipeline(query_fasta, target_fasta, mode = 'local', query_embedding_path = None, target_embedding_path = None, search_result_setting = None, output_path = None, if_stdout = True):
    print(f"resFIsher CPcheck2 {CPcheck2}\n")
    print(f"resFIsher CPcheck2 {CPcheck2}\n")
    print(f"resFIsher CPcheck2 {CPcheck2}\n")
    print(f"resFIsher CPcheck2 {CPcheck2}\n")
    print(f"resFIsher CPcheck2 {CPcheck2}\n")
    print(f"resFIsher CPcheck2 {CPcheck2}\n")

