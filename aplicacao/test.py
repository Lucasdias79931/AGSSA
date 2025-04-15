import os
import time
import csv
import platform
import psutil
import subprocess
import dill
import numpy as np
import shutil
from core_treino import pre_processamento_sequencias, executar_treinamento
from Bio import SeqIO

# Coletar informações da máquina
def get_cpu_name():
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if "model name" in line:
                    return line.split(":")[1].strip()
    except Exception:
        return platform.processor()

def get_motherboard_name():
    try:
        output = subprocess.check_output(["sudo", "dmidecode", "-t", "baseboard"], text=True)
        for line in output.splitlines():
            if "Product Name" in line:
                return line.split(":")[1].strip()
    except Exception:
        return f"Erro ao obter placa-mãe"

def get_memory_size():
    return f"{round(psutil.virtual_memory().total / (1024 ** 3), 2)} GB"

# Contar número de sequências
def analize_sequences_file(sequence_file_path):
    try:
        count = 0
        with open(sequence_file_path) as handle:
            for _, _ in SeqIO.FastaIO.SimpleFastaParser(handle):
                count += 1
        return count
    except Exception:
        return "Erro"

# Analisar modelo .obj com dill
def analisar_modelo(model_path):
    try:
        with open(model_path, "rb") as f:
            model = dill.load(f)
            return {
                "mean_resolution": round(float(np.sum(model.Resolution)), 5),
                "std_resolution": round(float(np.std(model.Resolution)), 5),
                "min_cluster_resolution": round(float(model.MinClusterResolution), 5),
                "n_clusters": model.NmbOfClusters,
                "n_classes": model.NmbOfClasses
            }
    except Exception as e:
        print(f"Erro ao carregar modelo: {e}")
        return {}

# Caminhos
here = os.path.abspath(os.path.dirname(__file__))
file_to_test_dir = os.path.join(here, "fileToTest")
results_dir = os.path.join(here, "resultTests")
agua_script = os.path.join(here, "AGUA", "AGUA_treinamento.py")
csv_path = os.path.join(results_dir, "resultados_teste.csv")

os.makedirs(results_dir, exist_ok=True)
especie = "spike"
data_teste = time.strftime("%d/%m/%Y %H:%M:%S")

# Info do sistema
dados_pc = {
    "cpu": get_cpu_name(),
    "memoria": get_memory_size(),
    "placa_mae": get_motherboard_name(),
    "data_execucao": data_teste
}

# Cabeçalho CSV
header = [
    "cpu", "memoria", "placa_mae", "data_execucao",
    "file_name", "especie", "numbers_of_sequence", "time_of_process",
    "mean_resolution", "std_resolution", "min_cluster_resolution",
    "n_clusters", "n_classes"
]

# Iniciar CSV
with open(csv_path, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)

    # Linha com informações do PC
    writer.writerow([
        dados_pc["cpu"], dados_pc["memoria"], dados_pc["placa_mae"], dados_pc["data_execucao"],
        "", "", "", "", "", "", "", "", ""
    ])

    # Início dos testes
    for i in range(1, 3):
        print(f"\n>> Testando conjunto {i}...")

        try:
            # Copiar arquivos para pasta de resultado
            original_seq = os.path.join(file_to_test_dir, f"_sequencias_treinamento.fasta")
            original_anot = os.path.join(file_to_test_dir, f"_anotacoes.txt")

            output_path = os.path.join(results_dir, f"result_{i}")
            os.makedirs(output_path, exist_ok=True)

            seq_path = os.path.join(output_path, "sequencias_treinamento.fasta")
            anot_path = os.path.join(output_path, "anotacoes.txt")

            shutil.copyfile(original_seq, seq_path)
            shutil.copyfile(original_anot, anot_path)

            start = time.time()
            sequencias_processadas = pre_processamento_sequencias(seq_path, anot_path, especie)
            executar_treinamento(agua_script, sequencias_processadas, anot_path, output_path, especie)
            tempo = time.time() - start

            model_path = os.path.join(output_path, f"model.{especie}.obj")
            metrics = analisar_modelo(model_path)

            writer.writerow([
                dados_pc["cpu"], dados_pc["memoria"], dados_pc["placa_mae"], dados_pc["data_execucao"],
                f"sequencias_spike{i}.fasta", especie,
                analize_sequences_file(original_seq),
                f"{tempo:.3f} s",
                metrics.get("mean_resolution", "N/A"),
                metrics.get("std_resolution", "N/A"),
                metrics.get("min_cluster_resolution", "N/A"),
                metrics.get("n_clusters", "N/A"),
                metrics.get("n_classes", "N/A")
            ])

        except Exception as e:
            print(f"Erro no teste {i}: {e}")
            writer.writerow([
                dados_pc["cpu"], dados_pc["memoria"], dados_pc["placa_mae"], dados_pc["data_execucao"],
                f"sequencias_spike{i}.fasta", especie,
                "Erro", "Erro", "Erro", "Erro", "Erro", "Erro", "Erro"
            ])
