import torch
cpu_num=5
torch.set_num_threads(cpu_num)
from plmcp.plmcp import plmCP_fisher_pipeline
plmCP_fisher_pipeline('example/rcsb_pdb_2PEL.fasta', 'example/rcsb_pdb_3CNA.fasta')
