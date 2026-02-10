#!/bin/bash
# docker_entrypoint.sh

echo "--- Setting up execution environment ---"

# 1. Copy User Files (INPUT)
# Pega TUDO que o usuário colocou na pasta e joga para dentro do sistema
if [ -d "/io" ]; then
    echo "Importing user files from 'input_and_output_files'..."
    # Copia para a raiz
    cp -f /io/* . 2>/dev/null || true
    # Copia TAMBÉM para a pasta arquivos (onde o executável procura)
    cp -f /io/* arquivos/ 2>/dev/null || true
else
    echo "WARNING: User directory not found. Using default internal examples."
fi

# 2. Run the Main Automation Script
echo "--- Starting Simulation (run_all.sh) ---"
./run_all.sh

# 3. Export Results (OUTPUT)
if [ -d "/io" ]; then
    echo "--- Exporting results to 'input_and_output_files' ---"
    
    # EM VEZ DE COPIAR SÓ 'saida.out', COPIAMOS TUDO QUE FOR RELEVANTE
    
    # Copia qualquer arquivo .out (saida.out, resultado.out, etc.)
    cp arquivos/*.out /io/ 2>/dev/null || true
    
    # Copia os gráficos
    cp arquivos/*.png /io/ 2>/dev/null || true
    
    # Copia logs e arquivos .txt gerados (caso ele mude a extensão)
    cp arquivos/mscdout.txt /io/ 2>/dev/null || true
    
    # Importante: Dá permissão total aos arquivos copiados para o usuário poder abrir/editar
    chmod 777 /io/* 2>/dev/null || true
    
    echo "Success! All results saved."
fi
