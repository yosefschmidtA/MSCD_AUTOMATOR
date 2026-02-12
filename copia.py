from flask import Flask, render_template, request, send_file, jsonify
import subprocess
import os
import zipfile
import io
from werkzeug.utils import secure_filename
import sys
import threading
import time

app = Flask(__name__)

# Variável global para saber se está rodando
simulacao_ativa = False

def salvar_arquivo_experimental(request):
    arquivo = request.files.get('experimentalFile')
    if arquivo and arquivo.filename != '':
        filename = secure_filename(arquivo.filename)
        arquivo.save(os.path.join(os.getcwd(), filename))
        arquivo.seek(0)
        arquivo.save(os.path.join("arquivos", filename))

@app.route('/')
def index():
    return render_template('index.html')

# --- FUNÇÃO QUE RODA EM SEGUNDO PLANO ---
def worker_simulacao():
    global simulacao_ativa
    log_file = "execution_log.txt"
    
    # Limpa o log antigo
    with open(log_file, "w") as f:
        f.write("--- Iniciando Simulação ---\n")

    # Comando com TEE: Mostra no seu terminal E salva no arquivo
    # 2>&1 garante que erros também vão para o log
    comando = f"bash run_all.sh 2>&1 | tee -a {log_file}"

    try:
        # shell=True permite usar o pipe '|' do tee
        subprocess.run(comando, shell=True, cwd=os.getcwd())
        
        # Marca no log que acabou
        with open(log_file, "a") as f:
            f.write("\n--- SIMULACAO CONCLUIDA ---")
            
    except Exception as e:
        with open(log_file, "a") as f:
            f.write(f"\nErro Crítico: {str(e)}")
            f.write("\n--- SIMULACAO CONCLUIDA (COM ERRO) ---")
    
    finally:
        simulacao_ativa = False

# --- ROTAS ---

@app.route('/iniciar_full', methods=['POST'])
def iniciar_full():
    global simulacao_ativa
    
    if simulacao_ativa:
        return "Já existe uma simulação rodando!", 409

    # 1. Preparação (Input e Arquivos)
    salvar_arquivo_experimental(request)
    conteudo_input = request.form.get('inputCluster')
    if not conteudo_input: return "Input vazio", 400
    
    with open("input_cluster.txt", "w") as f:
        f.write(conteudo_input)

    # 2. Inicia a Thread (Não trava o site)
    simulacao_ativa = True
    thread = threading.Thread(target=worker_simulacao)
    thread.daemon = True # Se fechar o server, a thread morre junto
    thread.start()

    return "Simulação iniciada", 200

@app.route('/stream_log')
def stream_log():
    """Lê o arquivo de log e retorna para o JS"""
    log_file = "execution_log.txt"
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            content = f.read()
        return content
    return "Aguardando log..."

@app.route('/check_status')
def check_status():
    """Informa ao JS se acabou"""
    return jsonify({"ativa": simulacao_ativa})

@app.route('/baixar_resultado_full')
def baixar_resultado_full():
    """Gera o ZIP apenas quando o usuário pedir no final"""
    memory_file = io.BytesIO()
    
    # Tenta achar o teory.out na raiz ou na pasta arquivos
    path_final = "teory.out"
    if not os.path.exists(path_final):
        path_final = os.path.join("arquivos", "teory.out")
    
    if os.path.exists(path_final):
        with zipfile.ZipFile(memory_file, 'w') as zf:
            zf.write(path_final, arcname="teory.out")
            if os.path.exists("execution_log.txt"):
                zf.write("execution_log.txt", arcname="log_execucao.txt")
        
        memory_file.seek(0)
        return send_file(memory_file, mimetype='application/zip', as_attachment=True, download_name='full_simulation_results.zip')
    else:
        return "Erro: Arquivo teory.out não encontrado. Verifique o log.", 404

# --- MANTIVE SUAS ROTAS ANTIGAS QUE FUNCIONAVAM BEM ---
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
    pass 

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
    pass

if __name__ == '__main__':
    app.run(debug=True, port=5000)
