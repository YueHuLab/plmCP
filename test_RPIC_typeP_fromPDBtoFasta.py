import torch
cpu_num=5
torch.set_num_threads(cpu_num)
from plmcp.plmcp import plmCP_fisher_pipeline
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1nk1a1.fasta', 'testExample/TransFastaByPDB/d1qdma1.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1d5ra1.fasta', 'testExample/TransFastaByPDB/d1rsy__.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1nls__.fasta', 'testExample/TransFastaByPDB/d2bqpa_.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1qasa2.fasta', 'testExample/TransFastaByPDB/d1rsy__.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1b6a_1.fasta', 'testExample/TransFastaByPDB/d1bia_1.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1b5ta_.fasta', 'testExample/TransFastaByPDB/d1k87a2.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1jwyb_.fasta', 'testExample/TransFastaByPDB/d1puja_.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1jwyb_.fasta', 'testExample/TransFastaByPDB/d1u0la2.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1nw5a_.fasta', 'testExample/TransFastaByPDB/d2adma_.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1gsa_1.fasta', 'testExample/TransFastaByPDB/d2hgsa1.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1qq5a_.fasta', 'testExample/TransFastaByPDB/d3chy__.fasta')
plmCP_fisher_pipeline('testExample/TransFastaByPDB/d1kiaa_.fasta', 'testExample/TransFastaByPDB/d1nw5a_.fasta')
