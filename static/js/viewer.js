/**
 * Visualizador 3D Robusto - MSCD Solver
 * Ignora o cabeçalho e foca apenas na estrutura após o primeiro ///
 */
function visualizarEstrutura() {
    const text = document.getElementById('inputCluster').value;
    const container = document.getElementById('container-3d');
    
    // 1. Localiza onde a estrutura realmente começa
    const firstDelimiterIndex = text.indexOf('///');
    
    if (firstDelimiterIndex === -1) {
        return alert("Erro: Nenhuma camada encontrada! Certifique-se de que a estrutura começa com '///'.");
    }

    // 2. Extrai apenas o texto estrutural (pula o cabeçalho)
    const structuralText = text.substring(firstDelimiterIndex + 3); 
    const layers = structuralText.split('///');

    // 3. Inicializa o visualizador
    let viewer = $3Dmol.createViewer(container, {backgroundColor: 'white'});
    viewer.clear();

    let currentZ = 0;
    const toRad = (deg) => deg * (Math.PI / 180);

    layers.forEach((layerContent) => {
        const lines = layerContent.trim().split('\n').map(l => l.trim());
        
        // Valida se o bloco de camada tem as informações necessárias
        if (lines.length < 9) return;

        const atomicNumber = parseInt(lines[0].split(' ')[0]);
        const aData = lines[2].split(/\s+/).map(Number);
        const bData = lines[4].split(/\s+/).map(Number);
        const oData = lines[6].split(/\s+/).map(Number);
        const spacing = parseFloat(lines[8]);

        // Acumula o Z para o empilhamento das camadas
        currentZ += spacing;

        // Conversão Polar -> Cartesiana
        const vecA = { x: aData[0] * Math.cos(toRad(aData[1])), y: aData[0] * Math.sin(toRad(aData[1])) };
        const vecB = { x: bData[0] * Math.cos(toRad(bData[1])), y: bData[0] * Math.sin(toRad(bData[1])) };
        const origin = { x: oData[0] * Math.cos(toRad(oData[1])), y: oData[0] * Math.sin(toRad(oData[1])) };

        // Cor da Prata (47)
        const color = atomicNumber === 47 ? '#C0C0C0' : '#4169E1';

        // Desenha a rede 3x3 para visualização do plano
        for (let i = -1; i <= 1; i++) {
            for (let j = -1; j <= 1; j++) {
                const x = origin.x + (i * vecA.x) + (j * vecB.x);
                const y = origin.y + (i * vecA.y) + (j * vecB.y);
                
                viewer.addSphere({
                    center: {x: x, y: y, z: currentZ},
                    radius: 0.6,
                    color: color
                });
            }
        }
    });

    viewer.zoomTo();
    viewer.render();
}
