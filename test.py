from Bio import SeqIO
import time
import sys
import os
import subprocess

def treinamento(filefasta, annotation, pythonFile, outputFile):
    try:
        command = ["python3", pythonFile, filefasta, annotation, outputFile, "spike"]
        subprocess.run(command, check=True)
    except IOError as e:
        print(f"IOError: {e}")
    except FileNotFoundError as e:
        print(f"Arquivo não encontrado: {e}")
    except Exception as e:
        print(f"Erro inesperado: {e}")

if len(sys.argv) < 3:
    print("Faltando parâmetros!")
    print("Uso: python3 treinamento.py <diretório com sequências> <número de arquivos>")
    sys.exit(1)

here = os.path.abspath(os.path.dirname(__file__))
resultDirectory = os.path.join(here, "resultTests")
os.makedirs(resultDirectory, exist_ok=True)

sequencias_dir = sys.argv[1]
num_arquivos = int(sys.argv[2])

for i in range(num_arquivos):
    filefasta = os.path.join(sequencias_dir, f"{i}.fasta")
    annotation = os.path.join(sequencias_dir, f"{i}.gff")  # ou outro caminho
    outputFile = os.path.join(resultDirectory, f"{i}_output.txt")
    pythonFile = "modelo.py"  # ou qualquer outro script que você deseje usar

    print(f"Treinando com {filefasta}...")
    treinamento(filefasta, annotation, pythonFile, outputFile)
