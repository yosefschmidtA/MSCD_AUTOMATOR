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
    """Roda o teo.py e move a imagem resultante para static (Com Debug Reforçado)"""
    try:
        # 1. Localiza o teory.out (Input)
        caminho_teory = "teory.out"
        if not os.path.exists(caminho_teory):
             # Tenta na pasta arquivos se não achar na raiz
             if os.path.exists("arquivos/teory.out"):
                caminho_teory = "arquivos/teory.out"
             else:
                return False, f"Arquivo teory.out não encontrado em: {os.getcwd()}"

        # 2. Prepara o comando
        # Importante: cwd=os.getcwd() garante que estamos na raiz /app do container
        cmd = ["python3", "teo.py", caminho_teory]
        
        # 3. Executa capturando TUDO
        resultado = subprocess.run(
            cmd,
            capture_output=True, # Pega o que o print() do python solta
            text=True,
            cwd=os.getcwd()
        )
        
        # Se o script quebrar (Exit Code != 0)
        if resultado.returncode != 0:
            return False, f"Crash no teo.py:\nSTDERR: {resultado.stderr}\nSTDOUT: {resultado.stdout}"

        # 4. Procura a imagem gerada (Output)
        nome_imagem = "grafico_polar3.png"
        
        # ONDE ESTÁ O WALLY? 
        # Vamos procurar a imagem na raiz E na pasta arquivos
        origem_encontrada = None
        
        if os.path.exists(nome_imagem):
            origem_encontrada = nome_imagem
        elif os.path.exists(os.path.join("arquivos", nome_imagem)):
            origem_encontrada = os.path.join("arquivos", nome_imagem)
            
        # 5. Move para a pasta Static (Pública)
        destino_final = "static/plot_resultado.png"

        if origem_encontrada:
            # Garante que a pasta static existe
            if not os.path.exists("static"):
                os.makedirs("static")
                
            # Move (substituindo se já existir)
            shutil.move(origem_encontrada, destino_final)
            
            return True, "Gráfico gerado com sucesso!"
        else:
            # AQUI ESTÁ O SEGREDO: Retorna o log para sabermos o que aconteceu
            msg_erro = (
                f"O script rodou (Exit 0), mas não gerou '{nome_imagem}'.\n"
                f"--- LOG DO SCRIPT ---\n{resultado.stdout}\n"
                f"--- ERROS SILENCIOSOS ---\n{resultado.stderr}"
            )
            return False, msg_erro

    except Exception as e:
        return False, f"Erro interno no servidor: {str(e)}"
