from flask import Blueprint, render_template, request, send_file, jsonify
import os
import io
import threading
import zipfile
from functools import wraps

# Importa nossos serviços
import file_service
import simulation_service
import iv_service

bp = Blueprint('main', __name__)

# O gunicorn roda com threads por causa do IV Curve. As rotas do MSCD gravam nos
# mesmos arquivos fixos (input_cluster.txt, arquivos/, os png de static/), e ate
# aqui so nao se atropelavam porque o worker sincrono atendia um pedido por vez.
# Esta trava mantem exatamente isso para elas; o IV Curve nao passa por ela.
_trava_mscd = threading.Lock()

def uma_por_vez(rota):
    @wraps(rota)
    def envolta(*args, **kwargs):
        with _trava_mscd:
            return rota(*args, **kwargs)
    return envolta

@bp.route('/')
def index():
    return render_template('home.html')

@bp.route('/mscd')
def mscd():
    return render_template('mscd.html')

@bp.route('/iv')
def iv():
    return render_template('iv.html')

# --- ROTAS IV CURVE (leitor de curvas IV de LEED) ---

def _erro_iv(erro, codigo=400):
    return jsonify({"erro": str(erro)}), codigo

@bp.route('/iv/upload', methods=['POST'])
def iv_upload():
    try:
        sessao = iv_service.cria_sessao(request.files.getlist('imagens'), request.files.getlist('xmls'))
        return jsonify(iv_service.info(sessao))
    except iv_service.ErroIV as erro:
        return _erro_iv(erro)

@bp.route('/iv/<sessao>/img/<int:energia>.png')
def iv_imagem(sessao, energia):
    try:
        return send_file(iv_service.png_cinza(sessao, energia), mimetype='image/png', max_age=3600)
    except iv_service.ErroIV as erro:
        return _erro_iv(erro, 404)

@bp.route('/iv/<sessao>/marca', methods=['POST'])
def iv_marca(sessao):
    d = request.get_json()
    try:
        return jsonify(iv_service.marca(sessao, int(d['energia']), float(d['x']), float(d['y'])))
    except iv_service.ErroIV as erro:
        return _erro_iv(erro)

@bp.route('/iv/<sessao>/calcula', methods=['POST'])
def iv_calcula(sessao):
    d = request.get_json()
    try:
        return jsonify(iv_service.calcula(sessao, d.get('ancoras'), d.get('metodo', 'fisico'),
                                          d.get('fundo', 'c'), int(d['largura']), d.get('centro')))
    except iv_service.ErroIV as erro:
        return _erro_iv(erro)

@bp.route('/iv/<sessao>/salva', methods=['POST'])
def iv_salva(sessao):
    d = request.get_json()
    try:
        nome, buf = iv_service.salva(sessao, d.get('nome'), d.get('ancoras'), d.get('metodo', 'fisico'),
                                     d.get('fundo', 'c'), int(d['largura']), bool(d.get('suavizar')),
                                     d.get('centro'))
    except iv_service.ErroIV as erro:
        return _erro_iv(erro)
    return send_file(buf, mimetype='text/plain', as_attachment=True, download_name=f'{nome}.txt')

# --- ROTAS FULL (ASSÍNCRONAS - TERMINAL WEB) ---

@bp.route('/rodar_full', methods=['POST'])
@uma_por_vez
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
@uma_por_vez
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
@uma_por_vez
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

@bp.route('/plotar_experimental', methods=['POST'])
@uma_por_vez
def plotar_experimental():
    # 1. Verifica se o navegador enviou o arquivo
    if 'file' not in request.files:
        return jsonify({"status": "error", "message": "Nenhum arquivo recebido pelo servidor."}), 400
    
    file = request.files['file']
    
    if file.filename == '':
        return jsonify({"status": "error", "message": "Nome do arquivo está vazio."}), 400

    # 2. Salva o arquivo temporariamente com um prefixo para NÃO brigar com o Fortran
    nome_temporario = "plot_temp_" + file.filename
    caminho_salvo = os.path.join(os.getcwd(), nome_temporario)
    file.save(caminho_salvo)

    # 3. Chama o serviço para rodar o exp.py neste arquivo temporário
    sucesso, msg = simulation_service.gerar_grafico_experimental(nome_temporario)
    
    # 4. Limpeza: apaga o arquivo temporário para não sujar o servidor
    if os.path.exists(caminho_salvo):
        os.remove(caminho_salvo)
    
    # 5. Retorna o resultado para o navegador
    if sucesso:
        import time
        return jsonify({
            "status": "success", 
            "url": f"/static/plot_exp.png?t={int(time.time())}"
        })
    else:
        return jsonify({"status": "error", "message": msg}), 500
