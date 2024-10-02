# plmCP

This is the implement of plmCP, a sequence-based, non-linear, flexible protein alignment tool for the detection of circular permutaion in proteins. 

## Install directly using Conda
conda env create -f environment.yml
  # activate
conda activate plmCP
## Install step by step using the Requirements file
Follow the steps in requirements_plmCP.sh

### Important 
Please give the directory of the ESM-1b model in the plmcp/embedding_generate.py file
Replace it.
# 
esm_model_path = '/home/data/t030413/.cache/torch/hub/checkpoints/esm1b_t33_650M_UR50S.pt'

##run tests
python test_site_2PEL_3CNA.py
python test_RPIC_typeP_fromPDBtoFasta.py
python test_site_RPIC.py
python test_site_RPIC_reverse.py

## Citation
Hu Y, Huang B. Zang C.Z. Detection of circular permutations by Protein Language Models.

##License
MIT license, See LICENSE file
The codes are based on the PLMalign (MIT license).Importantly, it change to the classical backtracking algorithm in Simth-Waterman algorithm.
##tips
The fasta file and fasta id file can be obtained from PDB file using the tools pdb2fasta.sh and pdb2Spit.sh in testExample/TransFastaAndID/ directory.
Please not that the fasta id file is not necessary needed. It was used to compare the result to the TMAlign and ICARUS conveniently with the animo acids ID.
In test_RPIC_typeP_fromPDBtoFasta.py, it give a example using only the fasta file.
The query and target can exchange their positions. For more details, please observe the different between test_site_RPIC.py and test_site_RPIC_reverse.py.
