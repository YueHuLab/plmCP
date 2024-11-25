import torch
cpu_num=5
torch.set_num_threads(cpu_num)
from plmcp.embedding_generate import esm_embedding_generate


def test_max_sequence_length_with_fasta(query_fasta):
    try:
        query_embeddings = esm_embedding_generate(query_fasta)
        print("Embedding generated successfully for the input fasta file.")
    except RuntimeError as e:
        print(f"Failed to generate embeddings. Error: {e}")

if __name__ == "__main__":
    # 假设 query_fasta 是你的输入 fasta 文件路径
    query_fasta = "5GJV_7fold.fasta"
    test_max_sequence_length_with_fasta(query_fasta)
