#!/bin/bash

# --- CONFIGURATION VARIABLES ---
# Defines the working folder name
PASTA_DE_TRABALHO="arquivos"
ARQUIVO_DE_ENTRADA="input_cluster.txt"

# --- START OF COMMAND SEQUENCE ---
echo "STARTING AUTOMATION..."
echo "-------------------------------------"

# 1. Checks if the input file exists in the main folder
if [ ! -f "$ARQUIVO_DE_ENTRADA" ]; then
    echo "ERROR: INPUT FILE '$ARQUIVO_DE_ENTRADA' NOT FOUND IN THE MAIN FOLDER."
    echo "SCRIPT WILL BE TERMINATED."
    exit 1
fi

# 2. Checks if the working folder exists.
if [ ! -d "$PASTA_DE_TRABALHO" ]; then
    echo "ERROR: DIRECTORY '$PASTA_DE_TRABALHO' NOT FOUND."
    echo "SCRIPT WILL BE TERMINATED."
    exit 1
fi

# 3. Copies the input file to the working folder
echo "COPYING FILE '$ARQUIVO_DE_ENTRADA' TO '$PASTA_DE_TRABALHO'..."
cp "$ARQUIVO_DE_ENTRADA" "$PASTA_DE_TRABALHO"/

# 4. Navigates to the working folder
cd "$PASTA_DE_TRABALHO"
echo "DIRECTORY CHANGED TO '$PASTA_DE_TRABALHO'."


# --- Bloco Novo (mexe no LD_LIBRARY_PATH E no PATH) ---
LIB_PATH=$(find / -name "libmpi_cxx.so.1" 2>/dev/null | head -n 1)

# Checks if the library was found
if [ -z "$LIB_PATH" ]; then
    echo "ERROR: Library 'libmpi_cxx.so.1' not found. Check your MPI installation."
    exit 1
fi

# Extrai o diretório da biblioteca (ex: /usr/lib64/mpi/gcc/openmpi/lib64)
LIB_DIR=$(dirname "$LIB_PATH")
# Deduz o diretório de binários (ex: /usr/lib64/mpi/gcc/openmpi/bin)
BIN_DIR=$(dirname "$LIB_DIR")/bin

# Exporta AMBOS os caminhos
export LD_LIBRARY_PATH="$LIB_DIR:$LD_LIBRARY_PATH"
export PATH="$BIN_DIR:$PATH"

echo "MPI LIBRARY DIR FOUND AT: $LIB_DIR"
echo "MPI BIN DIR ADDED TO PATH: $BIN_DIR"
make randmscd_parallel

# Verifica se o comando 'make' foi bem-sucedido
if [ $? -ne 0 ]; then
    echo "ERROR: 'make' command failed. Stopping script."
    exit 1
fi

echo "COMPILING FORTRAN EXECUTABLES..."


# --- Bloco de compilação Fortran ATUALIZADO ---
echo "COMPILING FORTRAN EXECUTABLES..."

# 1. Compilar phsh0 (sem alteração)
echo "Compiling phsh0..."
cd phsh0_src/
make FC=gfortran
cp phsh0 ../
cd ..

# 2. Compilar phsh1 (AGORA USANDO O MAKEFILE CORRETO)

# --- INÍCIO DO NOVO BLOCO 'poconv' (O CORRETO) ---
echo "COMPILING C++ EXECUTABLE 'poconv'..."

# 1. Entra no diretório fonte do poconv
cd poconv_src/

# 2. Roda a "fórmula mágica" que descobrimos
# (Adicionei -Wno-write-strings para limpar aqueles avisos de sintaxe antiga)
mpic++ -Wno-write-strings -o poconv poconv.cpp polation.cpp potentia.cpp userinfo.cpp userutil.cpp

# 3. Copia o executável final para a pasta principal
cp poconv ../

# 4. Volta para a pasta principal
cd ..

echo "COMPILATION COMPLETE."
# 5. Executes the cleanup of old files
echo "PERFORMING CLEANUP OF OLD FILES..."
find . -maxdepth 1 \( -name "ps*" -o -name "rm*" \) -not -name "rmconv" -not -name "rmconv.cpp" -not -name "rmconv.o" -not -name "psconv" -not -name "psconv.o" -not -name "psrm.x" -delete
rm output_header.txt
rm -f mufftin*.d
echo "CLEANUP COMPLETE."

# 6. Executes the Python script sequence
echo "EXECUTING PYTHON SCRIPTS..."
python3 gerador.py
sleep 3
python3 muff.py
sleep 3
python3 leitoF.py
sleep 3
python3 criador_final.py
sleep 3

# 7. Executes the randmscd program
echo "EXECUTING RANDMSCD PROGRAM..."
mpirun -np 5 randmscd_parallel output_header.txt

# 8. Copies the input file back to the main directory
echo "COPYING '$ARQUIVO_DE_ENTRADA' FILE BACK TO THE MAIN DIRECTORY..."
cp "$ARQUIVO_DE_ENTRADA" ../

# 9. Returns to the main directory
cd ..
echo "-------------------------------------"
echo "PROCESS COMPLETE. RETURNING TO THE MAIN DIRECTORY."
