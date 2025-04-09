from Bio import SeqIO
import time
import sys
import os
import subprocess

def treinamento(pythonFile, sequencesFile, annotation, outputFile, especie):
    try:
        command = ["python3", pythonFile, sequencesFile, annotation, outputFile, especie]
        subprocess.run(command, check=True)
    except IOError as e:
        print(f"IOError: {e}")
    except FileNotFoundError as e:
        print(f"Arquivo não encontrado: {e}")
    except Exception as e:
        print(f"Erro inesperado: {e}")

if len(sys.argv) < 4:
    print("Faltando parâmetros!")
    print("Uso: python3 treinamento.py <diretório com sequências> <número de arquivos>")
    
    print("Argumentos passados!")
    for arg in range(len(sys.argv)):
        print(f"{arg}: " + sys.argv[arg])
    sys.exit(1)

here = os.path.abspath(os.path.dirname(__file__))
resultDirectory = os.path.join(here, "resultTests")
os.makedirs(resultDirectory, exist_ok=True)

# executável de treinamento

pythonExc = sys.argv[1]

# Arquivo de sequencia

sequencias = sys.argv[2]

# Arquivo de anotação

anotacao = sys.argv[3]



especie = "spike"


treinamento(pythonExc, sequencias, anotacao, resultDirectory, especie)


