from Bio import SeqIO
import os
import sys
import hashlib

def transformInHash(sequence: str) -> str:
    return hashlib.sha512(sequence.encode()).hexdigest()

if __name__ == "__main__":
    # Caminho para as sequências
    try:
        arquivo_sequencias = sys.argv[1]
        print(arquivo_sequencias)

        pasta_upload = os.path.dirname(arquivo_sequencias)

        path_sequencias_filtradas = os.path.join(
            pasta_upload, 'sequencias_unicas.fasta')

        total_sequenciasProcessadas = 0
        hashes = set()

        with open(arquivo_sequencias, "r") as inputFile, open(path_sequencias_filtradas, "w") as outFile:
            for record in SeqIO.parse(inputFile, "fasta"):
                total_sequenciasProcessadas += 1
                
                seq = transformInHash(str(record.seq))

                if seq not in hashes:
                    hashes.add(seq)
                    SeqIO.write(record, outFile, "fasta")
    except FileNotFoundError:
        print(f"Erro: O arquivo '{arquivo_sequencias}' não foi encontrado.")
        sys.exit(1)
    except Exception as e:
        print(f"Erro inesperado: {e}")
        sys.exit(1)