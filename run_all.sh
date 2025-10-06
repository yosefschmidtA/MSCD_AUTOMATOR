#!/bin/bash

# --- VARIAVEIS DE CONFIGURACAO ---
# Define o nome da pasta de trabalho
PASTA_DE_TRABALHO="arquivos"
ARQUIVO_DE_ENTRADA="input_cluster.txt"

# --- INICIO DA SEQUENCIA DE COMANDOS ---
echo "INICIANDO A AUTOMACAO..."
echo "-------------------------------------"

# 1. Verifica se o arquivo de entrada existe na pasta principal
if [ ! -f "$ARQUIVO_DE_ENTRADA" ]; then
    echo "ERRO: O ARQUIVO DE ENTRADA '$ARQUIVO_DE_ENTRADA' NAO FOI ENCONTRADO NA PASTA PRINCIPAL."
    echo "O SCRIPT SERA ENCERRADO."
    exit 1
fi

# 2. Verifica se a pasta de trabalho existe e, se sim, navega para ela.
if [ ! -d "$PASTA_DE_TRABALHO" ]; then
    echo "ERRO: O DIRETORIO '$PASTA_DE_TRABALHO' NAO FOI ENCONTRADO."
    echo "O SCRIPT SERA ENCERRADO."
    exit 1
fi

# 3. Move o arquivo de entrada para a pasta de trabalho
echo "MOVENDO O ARQUIVO '$ARQUIVO_DE_ENTRADA' PARA '$PASTA_DE_TRABALHO'..."
cp "$ARQUIVO_DE_ENTRADA" "$PASTA_DE_TRABALHO"/

# 4. Navega para a pasta de trabalho
cd "$PASTA_DE_TRABALHO"
echo "DIRETORIO ALTERADO PARA '$PASTA_DE_TRABALHO'."

LIB_PATH=$(find / -name "libmpi_cxx.so.1" 2>/dev/null | head -n 1)

# Verifica se a biblioteca foi encontrada
if [ -z "$LIB_PATH" ]; then
    echo "ERRO: Biblioteca 'libmpi_cxx.so.1' não encontrada. Verifique a instalação do MPI."
    exit 1
fi

# Extrai o diretório do caminho completo e exporta
DIR_PATH=$(dirname "$LIB_PATH")
export LD_LIBRARY_PATH="$DIR_PATH:$LD_LIBRARY_PATH"

echo "BIBLIOTECA MPI ENCONTRADA EM: $DIR_PATH"
echo "VARIÁVEL LD_LIBRARY_PATH EXPORTADA."
# 5. Executa a limpeza de arquivos antigos
echo "REALIZANDO A LIMPEZA DE ARQUIVOS ANTIGOS..."
find . -maxdepth 1 \( -name "ps*" -o -name "rm*" \) -not -name "rmconv" -not -name "rmconv.cpp" -not -name "rmconv.o" -not -name "psconv" -not -name "psconv.o" -not -name "psrm.x" -delete
rm output_header.txt
rm -f mufftin*.d
echo "LIMPEZA CONCLUIDA."

# 6. Executa a sequência de scripts Python
echo "EXECUTANDO SCRIPTS PYTHON..."
python3 gerador.py
sleep 3
python3 muff.py
sleep 3
python3 leitoF.py
sleep 3
python3 criador_final.py
sleep 3

# 7. Executa o programa randmscd
echo "EXECUTANDO O PROGRAMA RANDMSCD..."
mpirun -np 10 randmscd_parallel output_header.txt

# 8. Copia o arquivo de entrada de volta para a pasta principal
echo "COPIANDO O ARQUIVO '$ARQUIVO_DE_ENTRADA' DE VOLTA PARA O DIRETORIO PRINCIPAL..."
cp "$ARQUIVO_DE_ENTRADA" ../

# 9. Retorna ao diretório principal
cd ..
echo "-------------------------------------"
echo "PROCESSO CONCLUIDO. RETORNANDO AO DIRETORIO PRINCIPAL."
