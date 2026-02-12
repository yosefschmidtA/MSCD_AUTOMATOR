from flask import Flask
from routes import bp  # Importa o Blueprint do arquivo routes.py

app = Flask(__name__)

# Registra as rotas
app.register_blueprint(bp)

if __name__ == '__main__':
    # Pode mudar para host='0.0.0.0' se for rodar no servidor para acesso externo
    app.run(debug=True, port=5000)
