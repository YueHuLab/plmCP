import torch
cpu_num=5
torch.set_num_threads(cpu_num)
from plmcp.plmcp import plmCP_fisher_pipelineWithAAidOnly
#plmCP_fisher_pipeline('test_1006/6u2m_C.fasta', 'test_1006/8b6s_A.fasta')
#plmCP_fisher_pipelineWithAAidOnly('/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/1r5m_A.id','/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/2d5l_A.id')
plmCP_fisher_pipelineWithAAidOnly('/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/2mta_H.id','/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/1kv9_A.id')
#plmCP_fisher_pipelineWithAAid('/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/1r5m_A.fasta','/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/1r5m_A.id','/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/2d5l_A.fasta','/home/data/t030413/foldseek/database/data_pdb/fastSinDup/readable/download2/2d5l_A.id')
