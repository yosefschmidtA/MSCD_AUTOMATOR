#!/bin/bash

# --- CONFIGURATION VARIABLES ---
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

# --- MPI Library Setup ---
LIB_PATH=$(find . -name "libmpi_cxx.so.1" 2>/dev/null | head -n 1)
if [ -z "$LIB_PATH" ]; then
    echo "ERROR: Library 'libmpi_cxx.so.1' not found. Check your MPI installation."
    exit 1
fi
DIR_PATH=$(dirname "$LIB_PATH")
export LD_LIBRARY_PATH="$DIR_PATH:$LD_LIBRARY_PATH"
echo "MPI LIBRARY FOUND AT: $DIR_PATH"
echo "LD_LIBRARY_PATH VARIABLE EXPORTED."

# --- START OF COMPILATION BLOCKS ---

# 1. Compilar randmscd_parallel (com limpeza manual)
echo "--- Cleaning and Compiling 'randmscd_parallel' ---"
# Remove o executável antigo e TODOS os arquivos .o do diretório principal
rm -f randmscd_parallel *.o
make randmscd_parallel
if [ $? -ne 0 ]; then
    echo "ERROR: 'make randmscd_parallel' command failed. Stopping script."
    exit 1
fi

# 2. Compilar phsh0 (com limpeza manual)
echo "--- Cleaning and Compiling 'phsh0' ---"
cd phsh0_src/
# Remove o executável antigo e TODOS os arquivos .o deste subdiretório
rm -f phsh0 *.o
make FC=gfortran
cp phsh0 ../
cd ..
if [ ! -f "phsh0" ]; then
    echo "ERROR: 'phsh0' compilation failed."
    exit 1
fi

# 3. Compilar poconv (com limpeza manual)
echo "--- Cleaning and Compiling 'poconv' ---"
cd poconv_src/
# Limpeza manual dos arquivos antigos
rm -f poconv *.o
# Compilação (esta linha é silenciosa se funcionar, mas mostrará erros se falhar)
mpic++ -Wno-write-strings -o poconv poconv.cpp polation.cpp potentia.cpp userinfo.cpp userutil.cpp
cp poconv ../
cd ..
if [ ! -f "poconv" ]; then
    echo "ERROR: 'poconv' compilation failed. Check for missing .h or .cpp files."
    exit 1
fi

echo "--- ALL COMPILATIONS COMPLETE ---"

# 5. Executes the cleanup of old *output* files
echo "PERFORMING CLEANUP OF OLD OUTPUT FILES..."
find . -maxdepth 1 \( -name "ps*" -o -name "rm*" \) \
 -not -name "psconv" -not -name "psconv.cpp" -not -name "psconv.o" \
 -not -name "psrm.x" -not -name "phsh1" -not -name "poconv" -delete
rm -f output_header.txt
rm -f mufftin*.d
echo "CLEANUP COMPLETE."

# 6. Executes the Python script sequence
echo "EXECUTING PYTHON SCRIPTS..."
python3 gerador.py
sleep 1
python3 muff.py
sleep 1
python3 leitoF.py
sleep 1
python3 criador_final.py
sleep 1

# 7. Executes the randmscd program
echo "EXECUTING RANDMSCD PROGRAM..."
mpirun -np 2 randmscd_parallel output_header.txt

# 8. Copies the input file back to the main directory
echo "COPYING '$ARQUIVO_DE_ENTRADA' FILE BACK TO THE MAIN DIRECTORY..."
cp "$ARQUIVO_DE_ENTRADA" ../

# 9. Returns to the main directory
cd ..
echo "-------------------------------------"
echo "PROCESS COMPLETE. RETURNING TO THE MAIN DIRECTORY."
