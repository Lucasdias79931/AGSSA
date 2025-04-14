import os
import subprocess
import datetime

def run_command(command):
    print(f"Executando: {command}")
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print("STDOUT:", stdout.decode())
    print("STDERR:", stderr.decode())
    return stdout.decode(), stderr.decode(), process.returncode


def pre_processamento_sequencias(sequencia_path, anotacoes_path, especie):
    aplicacao_path = os.path.abspath(os.path.dirname(__file__))
    referencia_path = os.path.join(aplicacao_path, 'database/sequencia_referencia')

    try:
        script1 = 'AGUA/pre-processamento/01_remover_sequencias_duplicadas.py'
        script2 = 'AGUA/pre-processamento/02_alinhar_sequencias.py'
        script3 = 'AGUA/pre-processamento/03_filtrar_sequencias_ruins.py'
        script4 = 'AGUA/pre-processamento/04_filtrar_anotacoes.py'

        referencia = os.path.join(referencia_path, f'sars-cov-2_spike.fasta')

        command1 = f'python3 {script1} {sequencia_path}'
        _, stderr1, code1 = run_command(command1)
        if code1 != 0:
            raise Exception(f"Erro no script 01: {stderr1}")

        output1 = sequencia_path.replace("sequencias_treinamento.fasta", "sequencias_unicas.fasta")

        alinhado_fasta = output1.replace("sequencias_unicas.fasta", "sequencias_alinhadas.fasta")
        command2 = f'python3 {script2} {output1} {referencia} {alinhado_fasta}'
        _, stderr2, code2 = run_command(command2)
        if code2 != 0:
            raise Exception(f"Erro no script 02: {stderr2}")

        command3 = f'python3 {script3} {alinhado_fasta}'
        _, stderr3, code3 = run_command(command3)
        if code3 != 0:
            raise Exception(f"Erro no script 03: {stderr3}")

        output3 = alinhado_fasta.replace("sequencias_alinhadas.fasta", "sequencias_tratadas.fasta")

        command4 = f'python3 {script4} {output3} {anotacoes_path}'
        _, stderr4, code4 = run_command(command4)
        if code4 != 0:
            raise Exception(f"Erro no script 04: {stderr4}")

        return output3  # ← ESSENCIAL

    except Exception as e:
        raise e



def executar_treinamento(agua_path, sequencia_path, anotacoes_path, result_path, especie, informacoes_path=""):
    try:
        command = f'python3 {agua_path} {sequencia_path} {anotacoes_path} {result_path} {especie}'
        os.system(command)

        if informacoes_path:
            with open(informacoes_path, 'a') as f:
                f.write(f'Data de Término: {datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")}\n')

        return "Feito"
    except Exception as e:
        raise e
