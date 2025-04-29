import os
import subprocess
import sys

def run_command(command):
    result = subprocess.run(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout.strip(), result.stderr.strip(), result.returncode

def ler_cabecalhos_fasta(caminho_fasta):
    with open(caminho_fasta) as f:
        return [linha.strip()[1:] for linha in f if linha.startswith(">")]

def restaurar_cabecalhos_em_memoria(cabecalhos_originais, caminho_saida):
    with open(caminho_saida) as f:
        linhas = f.readlines()

    novas_linhas = []
    index = 0

    for linha in linhas:
        if linha.startswith(">"):
            if index < len(cabecalhos_originais):
                novas_linhas.append(f">{cabecalhos_originais[index]}\n")
                index += 1
            else:
                novas_linhas.append(linha)
        else:
            novas_linhas.append(linha)

    with open(caminho_saida, "w") as f:
        f.writelines(novas_linhas)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Uso: python 02_alinhar_sequencias.py <arquivo_sequencias> <arquivo_referencia> <arquivo_saida>")
        sys.exit(1)

    arquivo_sequencias = sys.argv[1]
    arquivo_referencia = sys.argv[2]
    arquivo_saida = sys.argv[3]

    if not os.path.exists(arquivo_sequencias):
        print(f"Arquivo {arquivo_sequencias} não encontrado.")
        sys.exit(1)

    # Armazena os cabeçalhos em uma lista
    cabecalhos_originais = ler_cabecalhos_fasta(arquivo_sequencias)

    # Gera arquivo .sam temporário
    arquivo_sam = arquivo_sequencias.replace(".fasta", ".sam")
    command1 = f"minimap2 -a -x asm5 --cs --sam-hit-only --secondary=no -t 10 {arquivo_referencia} {arquivo_sequencias} -o {arquivo_sam}"
    stdout1, stderr1, returncode1 = run_command(command1)
    print(f"Comando 1: {command1}")
    print(f"stdout1: {stdout1}")
    print(f"stderr1: {stderr1}")

    if returncode1 != 0 or not os.path.exists(arquivo_sam):
        print(f"Erro ao gerar o arquivo .sam para {arquivo_sequencias}")
        print(f"Saída de erro: {stderr1}")
        sys.exit(1)

    # Executa o gofasta
    command2 = f"gofasta sam toMultiAlign -s {arquivo_sam} -t 10 --reference {arquivo_referencia} > {arquivo_saida}"
    stdout2, stderr2, returncode2 = run_command(command2)
    print(f"Comando 2: {command2}")
    print(f"stdout2: {stdout2}")
    print(f"stderr2: {stderr2}")

    if returncode2 == 0 and os.path.exists(arquivo_saida):
        print("Alinhamento concluído. Restaurando cabeçalhos originais...")
        restaurar_cabecalhos_em_memoria(cabecalhos_originais, arquivo_saida)
        print("Cabeçalhos restaurados com sucesso.")
    else:
        print(f"Erro ao gerar o arquivo de alinhamento {arquivo_saida}")
        print(f"Saída de erro: {stderr2}")
