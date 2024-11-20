import torch
cpu_num=5
torch.set_num_threads(cpu_num)
from plmcp.plmcp import plmCP_fisher_pipeline
plmCP_fisher_pipeline('test_1006/6u2m_C.fasta', 'test_1006/8b6s_A.fasta')
