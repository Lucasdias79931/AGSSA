import os
import time
import shutil
import dill
import numpy as np
from core_treino import pre_processamento_sequencias, executar_treinamento
import csv
import platform
import psutil
import subprocess

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
    except Exception as e:
        return f"Erro ao obter placa-mãe: {e}"

def get_memory_size():
    return f"{round(psutil.virtual_memory().total / (1024 ** 3), 2)} GB"

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

# Diretórios base
base_dir = os.path.abspath(os.path.dirname(__file__))
file_to_test_dir = os.path.join(base_dir, "fileToTest")
results_dir = os.path.join(base_dir, "resultTests")
os.makedirs(results_dir, exist_ok=True)

agua_script = os.path.join(base_dir, "AGUA", "AGUA_treinamento.py")
especie = "spike"
data_teste = time.strftime("%d/%m/%Y %H:%M:%S")

log_path = os.path.join(results_dir, "testes_resultados.txt")


dados_do_computador = {
    "cpu": get_cpu_name(),
    "memoria": get_memory_size(),
    "placa_mae": get_motherboard_name()
}

# Grava a introdução e informações do sistema uma única vez
with open(log_path, "a") as log:
    log.write("### TESTE AUTOMATIZADO AGSSA ###\n")
    log.write(f"Data de execução: {data_teste}\n")
    log.write("Informações da máquina:\n")
    for k, v in dados_do_computador.items():
        log.write(f"- {k}: {v}\n")
    log.write("\n")

# Loop para todos os pares de teste
for i in range(1, 3):
    print(f"\n>> Testando conjunto {i}...")

   
    original_seq = os.path.join(file_to_test_dir, f"_sequencias_treinamento.fasta")
    original_anot = os.path.join(file_to_test_dir, f"_anotacoes.txt")
    output_path = os.path.join(results_dir, f"result_{i}")
    os.makedirs(output_path, exist_ok=True)

    # Cópias para preservar os originais
    seq_path = os.path.join(output_path, "sequencias_treinamento.fasta")
    anot_path = os.path.join(output_path, "anotacoes.txt")
    shutil.copyfile(original_seq, seq_path)
    shutil.copyfile(original_anot, anot_path)

    info_path = os.path.join(output_path, "informacoes_processamento.txt")

    with open(info_path, "w") as f:
        f.write(f"Data de Início: {data_teste}\n")

    try:
        start = time.time()

        # Pré-processamento
        sequencias_processadas = pre_processamento_sequencias(seq_path, anot_path, especie)
        print(f"Saída do pré-processamento: {sequencias_processadas}")
        if not os.path.exists(sequencias_processadas):
            print("Arquivo de sequências tratadas NÃO encontrado.")
        elif os.path.getsize(sequencias_processadas) == 0:
            print("Arquivo de sequências tratadas está VAZIO.")

        # Treinamento (sem e-mail)
        executar_treinamento(
            agua_script,
            sequencias_processadas,
            anot_path,
            output_path,
            especie,
            informacoes_path=info_path
        )


        tempo = time.time() - start

        model_file = os.path.join(output_path, f"model.{especie}.obj")
        metrics = analisar_modelo(model_file)

        with open(log_path, "a") as log:
            log.write(f"\n--- TESTE {i} ---\n")
            log.write(f"Data: {data_teste}\n")
            log.write(f"Arquivo de entrada: sequencias_spike{i}.fasta\n")
            log.write(f"Tempo de execução: {tempo:.3f} segundos\n")
            for k, v in metrics.items():
                log.write(f"{k}: {v}\n")

    except Exception as e:
        print(f"Erro no teste {i}: {e}")
        with open(log_path, "a") as log:
            log.write(f"\n--- TESTE {i} ---\n")
            log.write(f"Erro: {e}\n")
