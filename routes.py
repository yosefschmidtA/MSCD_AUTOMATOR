from flask import Blueprint, render_template, request, send_file, jsonify
import os
import io
import zipfile

# Importa nossos serviços
import file_service
import simulation_service

bp = Blueprint('main', __name__)

@bp.route('/')
def index():
    return render_template('index.html')

# --- ROTAS FULL (ASSÍNCRONAS - TERMINAL WEB) ---

@bp.route('/rodar_full', methods=['POST'])
def iniciar_full():
    # 1. Salva inputs
    file_service.salvar_arquivo_experimental(request)
    conteudo = request.form.get('inputCluster')
    
    if not conteudo: 
        return "Input vazio", 400
    
    with open("input_cluster.txt", "w") as f:
        f.write(conteudo)

    # 2. Inicia o motor em segundo plano
    sucesso = simulation_service.iniciar_thread_full()
    
    if not sucesso:
        return "Já existe uma simulação rodando! Aguarde o término.", 409
        
    return "Iniciado", 200

@bp.route('/stream_log')
def stream_log():
    # Rota que o JS chama a cada 1 segundo
    return simulation_service.ler_log_atual()

@bp.route('/baixar_resultado_full')
def baixar_resultado_full():
    # Gera o ZIP apenas no final
    mem_file = file_service.gerar_zip_memoria("teory.out", "full_simulation.zip", log_file="execution_log.txt")
    if mem_file:
        return send_file(mem_file, mimetype='application/zip', as_attachment=True, download_name='full_simulation_results.zip')
    return "Erro: Arquivo teory.out não encontrado. Verifique o log no terminal.", 404

# --- ROTAS HALF (SÍNCRONAS - BOTÃO VERDE) ---
# (Mantenha igual ao que já estava funcionando)

@bp.route('/rodar', methods=['POST'])
def rodar_simulacao():
    file_service.salvar_arquivo_experimental(request)
    conteudo = request.form.get('inputCluster')
    
    if not conteudo: return "Input vazio", 400
    with open("input_cluster.txt", "w") as f: f.write(conteudo)

    sucesso, msg = simulation_service.rodar_half_script()
    if not sucesso: return f"Erro: {msg}", 500

    # Gera ZIP do Half
    memory_file = io.BytesIO()
    pasta_arquivos = "arquivos"
    with zipfile.ZipFile(memory_file, 'w') as zf:
        if os.path.exists(pasta_arquivos):
            for file in os.listdir(pasta_arquivos):
                if (file.startswith('ps') and file.count('.') >= 2) or \
                   (file.startswith('rm') and not file.startswith('psrm')):
                    zf.write(os.path.join(pasta_arquivos, file), file)
    
    memory_file.seek(0)
    return send_file(memory_file, mimetype='application/zip', as_attachment=True, download_name='resultados_mscd.zip')

@bp.route('/download_exemplo')
def download_exemplo():
    mem_file = file_service.gerar_zip_exemplos()
    return send_file(mem_file, mimetype='application/zip', as_attachment=True, download_name='exemplos_mscd.zip')

@bp.route('/gerar_grafico', methods=['POST'])
def gerar_grafico_rota():
    # Chama o serviço que roda o teo.py
    sucesso, msg = simulation_service.gerar_grafico()

    if sucesso:
        import time
        # Retorna o caminho da imagem com timestamp para não usar cache velho
        return jsonify({
            "status": "success",
            "url": f"/static/plot_resultado.png?t={int(time.time())}"
        })
    else:
        return jsonify({"status": "error", "message": msg}), 500
