import os
import subprocess
import sys

def run_command(command):
    result = subprocess.run(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout.strip(), result.stderr.strip(), result.returncode

if __name__ == "__main__":
    # Verifica se os argumentos necessários foram passados
    if len(sys.argv) < 4:
        print("Uso: python 02_alinhar_sequencias.py <arquivo_sequencias> <arquivo_referencia> <arquivo_saida>")
        sys.exit(1)

    # Argumentos esperados
    arquivo_sequencias = sys.argv[1]
    arquivo_referencia = sys.argv[2]
    arquivo_saida = sys.argv[3]

    if not os.path.exists(arquivo_sequencias):
        print(f"Arquivo {arquivo_sequencias} não encontrado.")
        sys.exit(1)

    # Gera nome para o .sam temporário
    arquivo_sam = arquivo_sequencias.replace(".fasta", ".sam")

    # Comando 1: minimap2
    command1 = f"minimap2 -a -x asm5 --cs --sam-hit-only --secondary=no -t 10 {arquivo_referencia} {arquivo_sequencias} -o {arquivo_sam}"
    stdout1, stderr1, returncode1 = run_command(command1)
    print(f"Comando 1: {command1}")
    print(f"stdout1: {stdout1}")
    print(f"stderr1: {stderr1}")

    if returncode1 != 0 or not os.path.exists(arquivo_sam):
        print(f"Erro ao gerar o arquivo .sam para {arquivo_sequencias}")
        print(f"Saída de erro: {stderr1}")
        sys.exit(1)

    # Comando 2: gofasta
    command2 = f"gofasta sam toMultiAlign -s {arquivo_sam} -t 10 --reference {arquivo_referencia} > {arquivo_saida}"
    stdout2, stderr2, returncode2 = run_command(command2)
    print(f"Comando 2: {command2}")
    print(f"stdout2: {stdout2}")
    print(f"stderr2: {stderr2}")

    if returncode2 == 0 and os.path.exists(arquivo_saida):
        print("Processamento concluído com sucesso.")
    else:
        print(f"Erro ao gerar o arquivo de alinhamento {arquivo_saida}")
        print(f"Saída de erro: {stderr2}")
