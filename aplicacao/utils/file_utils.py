import os


def ler_especies():
    especies_path = os.path.join(os.getcwd(), 'database/especies.txt')
    with open(especies_path, 'r') as file:
        especies = [line.strip() for line in file.readlines()]
    return especies


def conta_quantidade_sequencias(filepath):
    count_seqs = 0
    with open(filepath, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                count_seqs += 1
    return count_seqs
