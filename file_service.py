import os
import zipfile
import io
from werkzeug.utils import secure_filename

def salvar_arquivo_experimental(request):
    """Salva o arquivo experimental (se enviado) na raiz e na pasta arquivos"""
    arquivo = request.files.get('experimentalFile')
    if arquivo and arquivo.filename != '':
        filename = secure_filename(arquivo.filename)
        # Salva na raiz
        arquivo.save(os.path.join(os.getcwd(), filename))
        # Salva backup
        arquivo.seek(0)
        if not os.path.exists("arquivos"):
            os.makedirs("arquivos")
        arquivo.save(os.path.join("arquivos", filename))
        return True
    return False

def gerar_zip_memoria(arquivo_alvo, nome_zip, log_file=None):
    """Gera um ZIP na memória RAM pronto para download"""
    memory_file = io.BytesIO()
    
    # Define onde buscar o arquivo (raiz ou pasta arquivos)
    path_final = arquivo_alvo
    if not os.path.exists(path_final):
        path_final = os.path.join("arquivos", arquivo_alvo)
    
    if os.path.exists(path_final):
        with zipfile.ZipFile(memory_file, 'w') as zf:
            zf.write(path_final, arcname=arquivo_alvo)
            if log_file and os.path.exists(log_file):
                zf.write(log_file, arcname="log_execucao.txt")
        
        memory_file.seek(0)
        return memory_file
    return None

def gerar_zip_exemplos():
    """Gera o ZIP com os arquivos de exemplo"""
    memory_file = io.BytesIO()
    with zipfile.ZipFile(memory_file, 'w') as zf:
        for nome in ["inputexample.txt", "experimentalexample.txt", "exemplo.txt"]:
            caminho = os.path.join("arquivos", nome)
            if os.path.exists(caminho):
                zf.write(caminho, arcname=nome)
    memory_file.seek(0)
    return memory_file
