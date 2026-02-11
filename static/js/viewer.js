/**
 * Visualizador 3D Robusto - MSCD Solver
 * Features:
 * - Ignora cabeçalho
 * - Empilhamento Top-Down
 * - Cores CPK (1-80)
 * - Raio ajustado
 */

// Dicionário de Cores CPK (Elementos 1 ao 80)
const atomColors = {
    1: '#FFFFFF',  // H
    2: '#D9FFFF',  // He
    3: '#CC80FF',  // Li
    4: '#C2FF00',  // Be
    5: '#FFB5B5',  // B
    6: '#909090',  // C
    7: '#3050F8',  // N
    8: '#FF0D0D',  // O (Vermelho)
    9: '#90E050',  // F
    10: '#B3E2F5', // Ne
    11: '#AB5CF2', // Na
    12: '#8AFF00', // Mg
    13: '#BFA6A6', // Al
    14: '#F0C8A0', // Si
    15: '#FF8000', // P
    16: '#FFFF30', // S
    17: '#1FF01F', // Cl
    18: '#80D1E3', // Ar
    19: '#8F40D4', // K
    20: '#3DFF00', // Ca
    21: '#E6E6E6', // Sc
    22: '#BFC2C7', // Ti
    23: '#A6A6AB', // V
    24: '#8A99C7', // Cr
    25: '#9C7AC7', // Mn
    26: '#E06633', // Fe (Ferrugem)
    27: '#F090A0', // Co
    28: '#50D050', // Ni
    29: '#C88033', // Cu
    30: '#7D80B0', // Zn
    31: '#C28F8F', // Ga (Vinho/Vermelho claro - Importante para Ga2O3)
    32: '#668F8F', // Ge
    33: '#BD80E3', // As
    34: '#FFA100', // Se
    35: '#A62929', // Br
    36: '#5CB8D1', // Kr
    37: '#702EB0', // Rb
    38: '#00FF00', // Sr
    39: '#94FFFF', // Y
    40: '#94E0E0', // Zr
    41: '#73C2C9', // Nb
    42: '#54B5B5', // Mo
    43: '#3B9E9E', // Tc
    44: '#248F8F', // Ru
    45: '#0A7D8C', // Rh
    46: '#006985', // Pd
    47: '#C0C0C0', // Ag (Prata - Importante)
    48: '#FFD98F', // Cd
    49: '#A67573', // In
    50: '#668080', // Sn
    51: '#9E63B5', // Sb
    52: '#D47A00', // Te
    53: '#940094', // I
    54: '#429EB0', // Xe
    55: '#57178F', // Cs
    56: '#00C900', // Ba
    57: '#70D4FF', // La
    58: '#FFFFC7', // Ce
    59: '#D9FFC7', // Pr
    60: '#C7FFC7', // Nd
    61: '#A3FFC7', // Pm
    62: '#8FFFC7', // Sm
    63: '#61FFC7', // Eu
    64: '#45FFC7', // Gd
    65: '#30FFC7', // Tb
    66: '#1FFFC7', // Dy
    67: '#00FF9C', // Ho
    68: '#00E675', // Er
    69: '#00D452', // Tm
    70: '#00BF38', // Yb
    71: '#00AB24', // Lu
    72: '#4DC2FF', // Hf
    73: '#4DA6FF', // Ta
    74: '#2194D6', // W
    75: '#267DAB', // Re
    76: '#266696', // Os
    77: '#175487', // Ir
    78: '#D0D0E0', // Pt
    79: '#FFD123', // Au (Ouro)
    80: '#B8B8D0', // Hg
};

function visualizarEstrutura() {
    const text = document.getElementById('inputCluster').value;
    const container = document.getElementById('container-3d');
    
    // 1. Localiza onde a estrutura realmente começa
    const firstDelimiterIndex = text.indexOf('///');
    
    if (firstDelimiterIndex === -1) {
        return alert("Erro: Nenhuma camada encontrada! Certifique-se de que a estrutura começa com '///'.");
    }

    // 2. Extrai apenas o texto estrutural
    const structuralText = text.substring(firstDelimiterIndex + 3); 
    const layers = structuralText.split('///');

    // 3. Inicializa o visualizador
    let viewer = $3Dmol.createViewer(container, {backgroundColor: 'white'});
    viewer.clear();

    let currentZ = 0;
    const toRad = (deg) => deg * (Math.PI / 180);

    layers.forEach((layerContent) => {
        const lines = layerContent.trim().split('\n').map(l => l.trim());
        
        if (lines.length < 9) return;

        const atomicNumber = parseInt(lines[0].split(' ')[0]);
        const aData = lines[2].split(/\s+/).map(Number);
        const bData = lines[4].split(/\s+/).map(Number);
        const oData = lines[6].split(/\s+/).map(Number);
        const spacing = parseFloat(lines[8]);

        // Empilha para baixo (Layer 1 no topo)
        currentZ -= spacing;

        const vecA = { x: aData[0] * Math.cos(toRad(aData[1])), y: aData[0] * Math.sin(toRad(aData[1])) };
        const vecB = { x: bData[0] * Math.cos(toRad(bData[1])), y: bData[0] * Math.sin(toRad(bData[1])) };
        const origin = { x: oData[0] * Math.cos(toRad(oData[1])), y: oData[0] * Math.sin(toRad(oData[1])) };

        // Busca a cor no dicionário ou usa um magenta "de erro" se não achar
        const color = atomColors[atomicNumber] || '#FF00FF';

        for (let i = -2; i <= 2; i++) {
            for (let j = -2; j <= 2; j++) {
                const x = origin.x + (i * vecA.x) + (j * vecB.x);
                const y = origin.y + (i * vecA.y) + (j * vecB.y);
                
                viewer.addSphere({
                    center: {x: x, y: y, z: currentZ},
                    radius: 0.2, 
                    color: color
                });
            }
        }
    });

    viewer.zoomTo();
    viewer.render();
}
