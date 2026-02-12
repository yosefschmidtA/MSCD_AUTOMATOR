from flask import Flask, render_template, request, send_file
import subprocess
import os
import zipfile
import io

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/rodar', methods=['POST'])
def rodar_simulacao():
    conteudo_input = request.form.get('inputCluster')
    if not conteudo_input:
        return "Erro: O input está vazio!", 400

    # ESCREVE APENAS NO ARQUIVO DE ENTRADA NA RAIZ
    with open("input_cluster.txt", "w") as f:
        f.write(conteudo_input)

    try:
        print("--- Chamando seu script de execução (run_half.sh) ---")
        # O Python apenas inicia o seu script Bash. 
        # A inteligência de execução continua sendo o seu .sh original.
        subprocess.run(["bash", "run_half.sh"], check=True)

        # PREPARAÇÃO DO ZIP (Totalmente em memória RAM)
        memory_file = io.BytesIO()
        pasta_arquivos = "arquivos"
        
        with zipfile.ZipFile(memory_file, 'w') as zf:
            if os.path.exists(pasta_arquivos):
                for file in os.listdir(pasta_arquivos):
                    caminho_completo = os.path.join(pasta_arquivos, file)
                    
                    # FILTRO RÍGIDO DE SELEÇÃO (Apenas LEITURA):
                    # 1. Phase Shifts: ps1.1.txt, ps2.1.txt, etc. (Precisa ter 2 pontos)
                    is_ps_data = file.startswith('ps') and file.count('.') >= 2 and file.endswith('.txt')
                    
                    # 2. Radial Matrix: rm1.txt, rm2.txt (Ignora psrmin)
                    is_rm_data = file.startswith('rm') and file.endswith('.txt') and not file.startswith('psrm')

                    if is_ps_data or is_rm_data:
                        print(f"Lendo para o ZIP: {file}")
                        # O comando 'write' do zipfile apenas COPIA o conteúdo para o ZIP
                        zf.write(caminho_completo, file)

        memory_file.seek(0)
        
        return send_file(
            memory_file,
            mimetype='application/zip',
            as_attachment=True,
            download_name='resultados_mscd.zip'
        )

    except subprocess.CalledProcessError as e:
        return f"Erro durante a execução do script: {e}", 500
@app.route('/download_exemplo')
def download_exemplo():
    # Cria um arquivo ZIP na memória (RAM)
    memory_file = io.BytesIO()

    # Abre o ZIP para escrita
    with zipfile.ZipFile(memory_file, 'w') as zf:
        # Adiciona o exemplo.txt
        caminho_1 = os.path.join("arquivos", "inputexample.txt")
        if os.path.exists(caminho_1):
            zf.write(caminho_1, arcname="inputexample.txt")

        # Adiciona o experimentalexample.txt
        caminho_2 = os.path.join("arquivos", "experimentalexample.txt")
        if os.path.exists(caminho_2):
            zf.write(caminho_2, arcname="experimentalexample.txt")

    # Volta o ponteiro para o início do arquivo na memória
    memory_file.seek(0)

    return send_file(
        memory_file,
        mimetype='application/zip',
        as_attachment=True,
        download_name='exemplos_mscd.zip'
    )
if __name__ == '__main__':
    # Roda em localhost. O debug=True ajuda a ver erros no terminal.
    app.run(debug=True, port=5000)
