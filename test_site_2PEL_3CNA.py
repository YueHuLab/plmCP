import torch
cpu_num=5
torch.set_num_threads(cpu_num)
from plmcp.plmcp import plmCP_fisher_pipelineWithAAid
plmCP_fisher_pipelineWithAAid('testExample/TransFastaAndID/test/2pel.fasta', 'testExample/TransFastaAndID/test/3cna.fasta','testExample/TransFastaAndID/test/2pel.id', 'testExample/TransFastaAndID/test/3cna.id')
