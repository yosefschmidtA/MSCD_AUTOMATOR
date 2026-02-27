import subprocess
import os
import threading
import shutil
import signal

# Configurações Globais
LOG_FILE = "execution_log.txt"
simulacao_ativa = False

def worker_simulacao():
    """Função que roda em segundo plano (Thread) com Monitoramento em Tempo Real"""
    global simulacao_ativa
    
    # 1. Limpa/Cria o log vazio com cabeçalho
    with open(LOG_FILE, "w") as f:
        f.write("--- Initiating Simulation (Waiting for processing) ---\n")

    # 2. O Comando (SEM O TEE)
    # Removemos o '| tee' porque o Python vai ler a saída e escrever no arquivo ele mesmo.
    # '2>&1' garante que erros também venham para a leitura.
    comando = "bash run_all.sh 2>&1"

    processo = None

    try:
        # Inicia o processo permitindo leitura em tempo real (Popen)
        processo = subprocess.Popen(
            comando, 
            shell=True, 
            cwd=os.getcwd(), 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True, 
            bufsize=1 # Importante: Buffer linha a linha
        )
        
        erro_detectado = False

        # Abre o arquivo de log para ir gravando o que acontece
        with open(LOG_FILE, "a") as f:
            # Loop que lê cada linha que o run_all.sh cospe
            for linha in processo.stdout:
                f.write(linha)
                f.flush() # Força gravar no disco pro JS ler rápido
                
                # --- O GUARDA DE TRÂNSITO ---
                if "Error 201" in linha:
                    f.write("\n\n⛔ DETECTED ERROR 201. KILLING PROCESS IMMEDIATELY.\n")
                    f.write("--- SIMULATION ABORTED (FILE NOT FOUND) ---\n")
                    erro_detectado = True
                    
                    # Mata o processo do run_all.sh
                    processo.kill() 
                    break # Sai do loop
        
        # Se saiu do loop, espera limpar o processo
        if processo:
            processo.wait()

        # Só escreve "COMPLETE" se NÃO detectou o erro 201
        if not erro_detectado:
            with open(LOG_FILE, "a") as f:
                f.write("\n--- Simulation Complete ---")
            
    except Exception as e:
        with open(LOG_FILE, "a") as f:
            f.write(f"\nCritical Python Error: {str(e)}")
            f.write("\n--- SIMULATION COMPLETE (WITH ERROR) ---")
    
    finally:
        # Libera a variável para permitir nova simulação
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
    return "Initializing log..."

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
def gerar_grafico_experimental(nome_arquivo_dinamico):
    """Roda o exp.py isolado do Fortran para não interferir na simulação"""
    try:
        caminho_exp = os.path.join(os.getcwd(), nome_arquivo_dinamico)
        
        if not os.path.exists(caminho_exp):
            return False, f"Arquivo {nome_arquivo_dinamico} não encontrado para plotar."

        # Roda o exp.py passando o nome temporário
        cmd = ["python3", "exp.py", caminho_exp]
        resultado = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd())

        # Se o script falhar
        if resultado.returncode != 0:
            return False, f"Crash no exp.py:\n{resultado.stderr}\nSTDOUT: {resultado.stdout}"

        # Procura a imagem gerada
        nome_imagem = "grafico_exp.png" 
        destino_final = "static/plot_exp.png"

        if os.path.exists(nome_imagem):
            if not os.path.exists("static"): 
                os.makedirs("static")
            
            # Move e sobrescreve a imagem antiga em static/
            shutil.move(nome_imagem, destino_final)
            return True, "Gráfico experimental gerado com sucesso!"
        else:
            return False, f"Script rodou, mas não gerou '{nome_imagem}'. Log: {resultado.stdout}"

    except Exception as e:
        return False, f"Erro interno: {str(e)}"
