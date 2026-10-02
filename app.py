from flask import Flask
from routes import bp  # Importa o Blueprint do arquivo routes.py

app = Flask(__name__)
# O leitor de IV recebe a pasta inteira de uma vez: 250 energias com .jpg e .xml
# ja passam das 1000 partes que o Werkzeug aceita por padrao.
app.config['MAX_FORM_PARTS'] = 20000

# Registra as rotas
app.register_blueprint(bp)

if __name__ == '__main__':
    # Pode mudar para host='0.0.0.0' se for rodar no servidor para acesso externo
    app.run(debug=True, port=5000)
