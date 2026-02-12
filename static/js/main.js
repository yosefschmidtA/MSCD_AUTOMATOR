/**
 * main.js - Gerencia Simulação e Terminal Ao Vivo
 */

let logInterval = null;

function executarSimulacao(rota) {
    const content = document.getElementById('inputCluster').value;
    const fileInput = document.getElementById('experimentalFile');
    const status = document.getElementById('status');
    const resultsCard = document.getElementById('resultsCard'); // Card do Run Half
    const downloadLink = document.getElementById('downloadLink'); // Link do Run Half

    // Validação
    if(!content.trim()) return alert("Erro: O campo de input está vazio!");

    // Prepara FormData
    const formData = new FormData();
    formData.append('inputCluster', content);
    if (fileInput && fileInput.files.length > 0) {
        formData.append('experimentalFile', fileInput.files[0]);
    }

    // === CENÁRIO 1: RUN HALF (Botão Verde) ===
    if (rota === '/rodar') {
        // UI Setup
        document.getElementById('terminal-container').style.display = 'none';
        document.getElementById('download-area-full').style.display = 'none';
        resultsCard.style.display = 'none';
        status.style.display = 'block';
        status.innerText = "⚙️ Gerando ps & rm...";

        fetch(rota, { method: 'POST', body: formData })
        .then(async response => {
            if (!response.ok) throw new Error(await response.text());
            return response.blob();
        })
        .then(blob => {
            const url = window.URL.createObjectURL(blob);
            downloadLink.href = url;
            downloadLink.download = 'mscd_results.zip';
            
            status.style.display = 'none';
            resultsCard.style.display = 'block';
        })
        .catch(err => {
            alert("Erro: " + err.message);
            status.style.display = 'none';
        });
        return;
    }

    // === CENÁRIO 2: RUN FULL (Botão Laranja) ===
    if (rota === '/rodar_full') {
        // UI Setup
        resultsCard.style.display = 'none';
        status.style.display = 'none';
        
        const terminalContainer = document.getElementById('terminal-container');
        const terminalLog = document.getElementById('terminal-log');
        const downloadArea = document.getElementById('download-area-full');

        terminalContainer.style.display = 'block';
        downloadArea.style.display = 'none';
        terminalLog.innerText = "🚀 Iniciando processo...\n";

        // 1. Dispara o início da simulação
        fetch('/rodar_full', { method: 'POST', body: formData })
        .then(async response => {
            if (!response.ok) {
                // Se já tiver rodando (409) ou erro (500)
                const msg = await response.text();
                throw new Error(msg);
            }
            // 2. Se deu 200 OK, começa a ler o log
            iniciarLeituraDoLog();
        })
        .catch(err => {
            alert("Não foi possível iniciar: " + err.message);
            terminalLog.innerText += "\n[ERRO AO INICIAR]: " + err.message;
        });
    }
}

function iniciarLeituraDoLog() {
    const terminalLog = document.getElementById('terminal-log');
    
    // Garante que não temos 2 intervalos rodando ao mesmo tempo
    if (logInterval) clearInterval(logInterval);

    logInterval = setInterval(() => {
        fetch('/stream_log')
        .then(res => res.text())
        .then(textoLog => {
            // Atualiza o texto da tela preta
            terminalLog.innerText = textoLog;
            
            // Rola para o final (auto-scroll)
            terminalLog.scrollTop = terminalLog.scrollHeight;

            // Verifica se acabou
            if (textoLog.includes("--- SIMULACAO CONCLUIDA ---")) {
                clearInterval(logInterval); // Para de ler
                document.getElementById('download-area-full').style.display = 'block'; // Mostra botão de baixar
            }
        })
        .catch(err => console.error("Erro ao ler log:", err));
    }, 100); // Lê a cada 1 segundo (1000ms)
}

function plotarGrafico() {
    const status = document.getElementById('status');
    const plotArea = document.getElementById('plot-area');
    const plotImg = document.getElementById('plot-img');

    // Mostra que está carregando (reaproveitando a div de status ou criando um alerta)
    if(status) {
        status.style.display = 'block';
        status.innerText = "🎨 Gerando gráfico...";
    }

    // Esconde a área antiga enquanto carrega a nova
    if(plotArea) plotArea.style.display = 'none';

    fetch('/gerar_grafico', { method: 'POST' })
    .then(res => res.json())
    .then(data => {
        if(status) status.style.display = 'none';

        if (data.status === 'success') {
            if(plotImg && plotArea) {
                // O timestamp (?t=...) força o navegador a baixar a imagem nova
                plotImg.src = data.url;
                plotArea.style.display = 'block';
            }
        } else {
            alert("Erro ao plotar: " + data.message);
        }
    })
    .catch(err => {
        if(status) status.style.display = 'none';
        alert("Erro de comunicação ao tentar plotar: " + err);
    });
}
