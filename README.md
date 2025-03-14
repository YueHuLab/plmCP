# plmCP

This is the implement of plmCP, a sequence-based, non-linear, flexible protein alignment tool for the detection of circular permutaion in proteins. 

## Install directly using Conda
conda env create -f environment.yml
  # activate
conda activate plmCP
## Install step by step using the Requirements file
Follow the steps in requirements_plmCP.sh

### Important 
Please give the directory of the ESM-1b, ESM-2, ProtT5 model in the plmcp/embedding_generate.py file
Replace it. There are four branches now.
# esm1b (esm_embedding_generate function)
esm_model_path = '/home/data/t030413/.cache/torch/hub/checkpoints/esm1b_t33_650M_UR50S.pt'
# https://github.com/YueHuLab/plmCP/tree/main
# or esm2 (esm_embedding_generate function)
esm_model_path = '/home/data/t030413/.cache/torch/hub/checkpoints/esm2_t33_650M_UR50D.pt'
#  https://github.com/YueHuLab/plmCP/tree/esm2
# or ProtT5 (prottrans_embedding_generate fucntion) and for ProtT5 you should also change the embedding function in plmCP.py
prottrans_model_path='/home/data/t030413/PLMalign/data/prot_t5_xl_uniref50' 
# https://github.com/YueHuLab/plmCP/tree/ProtT5
# we also build a version to extend the length of tokens by slide window, inspired by trRosetta
# https://github.com/YueHuLab/plmCP/tree/SlideWindow
# online running by Code Ocean
# https://codeocean.com/capsule/8856481/tree/v1

## run tests
python test_site_2PEL_3CNA.py
#
python test_RPIC_typeP_fromPDBtoFasta.py
#
python test_site_RPIC.py
#
python test_site_RPIC_reverse.py
## slide windows hyperparameter test
We tested different window lengths (300, 400 and 500) and shift lengths (100, 200 and 300) and found that the alignment length and score were generally consistent, demonstrating that our method has robustness within a certain range.This test was conducted using 2mta_H and 1kv9_A, and can be run with test_1201_slide.py. The length of 1kv9_A is 664, and after duplication, it reaches 1328, which will trigger the slide windows setting. This program uses a double-window strategy to slide over different parts of the sequence. By combining each pair of windows, it extracts features to capture long-range interactions within the sequence. Each position is accessed multiple times to enhance the robustness of feature representation. The final result is obtained by averaging the features from multiple windows to achieve a more stable representation. The double-window approach allows the model to capture interactions between different parts of the sequence, enhancing the representation by providing richer context information.



### most simple test and no residue id information
python test_2PEL_3CNA.py

## Citation
Hu, Yue, Bin Huang, Chun Zi Zang, and Jia Jie Xu. "Detection of circular permutations by Protein Language Models." Computational and Structural Biotechnology Journal 27 (2025): 214-220.

## License
MIT license, See LICENSE file
The codes are based on the PLMalign (MIT license).Importantly, it change to the classical backtracking algorithm in Simth-Waterman algorithm.
## tips
The fasta file and fasta id file can be obtained from PDB file using the tools pdb2fasta.sh and pdb2Spit.sh in testExample/TransFastaAndID/ directory.
Please not that the fasta id file is not necessary needed. It was used to compare the result to the TMAlign and ICARUS conveniently with the animo acids ID.
In test_RPIC_typeP_fromPDBtoFasta.py, it give a example using only the fasta file.
The query and target can exchange their positions. For more details, please observe the different between test_site_RPIC.py and test_site_RPIC_reverse.py.
## run it online on CodeOcean
We provide a code ocean version to run it by yourself.
https://codeocean.com/capsule/8856481/tree/v1
