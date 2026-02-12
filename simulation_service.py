import subprocess
import os
import threading
import shutil

# Configurações Globais
LOG_FILE = "execution_log.txt"
simulacao_ativa = False

def worker_simulacao():
    """Função que roda em segundo plano (Thread)"""
    global simulacao_ativa
    
    # 1. Limpa/Cria o log vazio com cabeçalho
    with open(LOG_FILE, "w") as f:
        f.write("--- Iniciando Simulação (Aguarde o processamento) ---\n")

    # 2. Comando Mágico:
    # "2>&1" pega erros; "| tee -a" mostra no terminal E salva no arquivo
    comando = f"bash run_all.sh 2>&1 | tee -a {LOG_FILE}"

    try:
        # shell=True é necessário para usar o pipe '|'
        subprocess.run(comando, shell=True, cwd=os.getcwd())
        
        # Marca o fim no log para o JS saber que acabou
        with open(LOG_FILE, "a") as f:
            f.write("\n--- SIMULACAO CONCLUIDA ---")
            
    except Exception as e:
        with open(LOG_FILE, "a") as f:
            f.write(f"\nErro Crítico no Python: {str(e)}")
            f.write("\n--- SIMULACAO CONCLUIDA (COM ERRO) ---")
    
    finally:
        simulacao_ativa = False

def iniciar_thread_full():
    """Tenta iniciar a thread. Retorna True se sucesso."""
    global simulacao_ativa
    
    if simulacao_ativa:
        return False
    
    simulacao_ativa = True
    thread = threading.Thread(target=worker_simulacao)
    thread.daemon = True # Se fechar o servidor, a thread morre junto
    thread.start()
    return True

def ler_log_atual():
    """Lê o conteúdo atual do arquivo de log"""
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, 'r') as f:
            return f.read()
    return "Inicializando log..."

def check_ativa():
    return simulacao_ativa

# Função Síncrona para o botão verde (Half)
def rodar_half_script():
    try:
        subprocess.run(["bash", "run_half.sh"], check=True)
        return True, "Sucesso"
    except subprocess.CalledProcessError as e:
        return False, str(e)
def gerar_grafico():
    """Roda o teo.py e move a imagem resultante para static"""
    try:
        # 1. Verifica se o teory.out existe (na raiz ou em arquivos)
        caminho_teory = "teory.out"
        if not os.path.exists(caminho_teory):
             caminho_teory = "arquivos/teory.out"
        
        if not os.path.exists(caminho_teory):
            return False, "Arquivo teory.out não encontrado para plotagem."

        # 2. Roda o script teo.py
        # Ele vai gerar o arquivo 'grafico_polar3.png' na raiz
        cmd = ["python3", "teo.py", caminho_teory]
        
        resultado = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )
        
        if resultado.returncode != 0:
            return False, f"Erro no script teo.py: {resultado.stderr}"

        # 3. Move a imagem gerada para a pasta public (static)
        # O nome original no seu script é 'grafico_polar3.png'
        nome_original = "grafico_polar3.png"
        destino_final = "static/plot_resultado.png"

        if os.path.exists(nome_original):
            # Garante que a pasta static existe
            if not os.path.exists("static"):
                os.makedirs("static")
                
            # Move e renomeia
            shutil.move(nome_original, destino_final)
            return True, "Gráfico gerado com sucesso!"
        else:
            return False, f"O script rodou, mas não gerou o arquivo '{nome_original}'."

    except Exception as e:
        return False, f"Erro interno: {str(e)}"
