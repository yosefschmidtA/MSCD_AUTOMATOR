// IV Curve: a janela do leed_iv.py no navegador. O calculo fica no servidor
// (leed_iv_core.py); aqui so a tela, os cliques e o estado do spot atual.

const RE_IMG = /^\d+\.(jpe?g|tiff?)$/i;
const RE_XML = /^\d+\.xml$/i;
// Cores categoricas fixas dos spots ja salvos (paleta validada do skill de dataviz).
const CORES_SALVOS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300'];
const COR_ATUAL = '#1f2328';
const COR_ENERGIA = '#d03b3b';
const COR_ANCORA = '#b8860b';

const $ = id => document.getElementById(id);

const S = {
    sessao: null, energias: [], largura: 15, forma: [0, 0], xml: {}, temCorrente: false,
    ancoras: {}, metodo: 'fisico', fundo: 'c', normalizar: false, suavizar: false,
    centro: null, centroAjustado: null, pos: {}, inten: {}, salvos: [], mensagem: '', i: 0,
    pedido: 0,
};
const cacheImg = new Map();  // energia -> Promise<{w, h, cinza, tela}>

// ------------------------------------------------------------------ envio das pastas

function arquivosValidos(input, re) {
    return Array.from(input.files).filter(f => re.test(f.name));
}

function atualizaContagem() {
    const imgs = arquivosValidos($('iv-pasta-img'), RE_IMG);
    const xmls = arquivosValidos($('iv-pasta-xml'), RE_XML).concat(arquivosValidos($('iv-pasta-img'), RE_XML));
    $('iv-conta-img').textContent = $('iv-pasta-img').files.length
        ? `${imgs.length} images found` : 'No folder selected';
    $('iv-conta-xml').textContent = $('iv-pasta-xml').files.length || xmls.length
        ? `${xmls.length} XML files found` : 'No folder selected';
    $('iv-carrega').disabled = imgs.length === 0;
}

function status(texto, erro = false) {
    $('iv-status').textContent = texto;
    $('iv-status').classList.toggle('iv-erro', erro);
}

function carrega() {
    const fd = new FormData();
    arquivosValidos($('iv-pasta-img'), RE_IMG).forEach(f => fd.append('imagens', f, f.name));
    arquivosValidos($('iv-pasta-img'), RE_XML).forEach(f => fd.append('xmls', f, f.name));
    // Na pasta propria, o XML vence o que veio junto das imagens: vai por ultimo.
    arquivosValidos($('iv-pasta-xml'), RE_XML).forEach(f => fd.append('xmls', f, f.name));

    const xhr = new XMLHttpRequest();
    xhr.open('POST', '/iv/upload');
    xhr.responseType = 'json';
    $('iv-carrega').disabled = true;
    $('iv-progresso').hidden = false;
    status('Uploading...');
    xhr.upload.onprogress = ev => {
        if (ev.lengthComputable) $('iv-progresso-barra').style.width = `${100 * ev.loaded / ev.total}%`;
    };
    xhr.onload = () => {
        $('iv-carrega').disabled = false;
        $('iv-progresso').hidden = true;
        $('iv-progresso-barra').style.width = '0';
        if (xhr.status !== 200) {
            status((xhr.response && xhr.response.erro) || `Upload failed (HTTP ${xhr.status}).`, true);
            return;
        }
        inicia(xhr.response);
    };
    xhr.onerror = () => {
        $('iv-carrega').disabled = false;
        $('iv-progresso').hidden = true;
        status('Upload failed: no response from the server.', true);
    };
    xhr.send(fd);
}

function inicia(r) {
    cacheImg.clear();
    Object.assign(S, {
        sessao: r.sessao, energias: r.energias, largura: r.largura, forma: r.forma, xml: r.xml,
        temCorrente: r.tem_corrente, ancoras: {}, normalizar: r.tem_corrente, centro: null,
        centroAjustado: null, pos: {}, inten: {}, salvos: [], i: 0,
        mensagem: 'Click on the center of the spot you want to measure.',
    });
    let txt = `${r.energias.length} energies loaded (${r.energias[0]} to ${r.energias.at(-1)} eV).`;
    if (r.sem_xml.length) {
        txt += ` ${r.sem_xml.length} without beam current in the XML: normalization by current disabled.`;
    }
    status(txt);
    $('iv-slider').max = r.energias.length - 1;
    $('iv-slider').value = 0;
    $('iv-ck-corrente').checked = S.normalizar;
    $('iv-ck-corrente-label').hidden = !S.temCorrente;
    $('iv-ferramenta').hidden = false;
    vaiPara(0);
    // Pre-carrega as outras para o passeio pelas energias nao travar.
    S.energias.forEach(e => imagem(e));
}

// ------------------------------------------------------------------ imagens

function percentil(hist, n, p) {
    const alvo = p / 100 * (n - 1);
    let acum = 0;
    for (let v = 0; v < 256; v++) {
        acum += hist[v];
        if (acum > alvo) return v;
    }
    return 255;
}

function imagem(e) {
    if (!cacheImg.has(e)) {
        cacheImg.set(e, new Promise((ok, falha) => {
            const im = new Image();
            im.onload = () => {
                const w = im.naturalWidth, h = im.naturalHeight;
                const cv = document.createElement('canvas');
                cv.width = w; cv.height = h;
                const cx = cv.getContext('2d');
                cx.drawImage(im, 0, 0);
                const dados = cx.getImageData(0, 0, w, h);
                const cinza = new Uint8Array(w * h);
                const hist = new Uint32Array(256);
                for (let k = 0; k < w * h; k++) {
                    cinza[k] = dados.data[4 * k];
                    hist[cinza[k]]++;
                }
                // Contraste so da tela, como o tvscl: nao mexe nos numeros da curva.
                const lo = percentil(hist, w * h, 1);
                const hi = Math.max(percentil(hist, w * h, 99.8), 1, lo + 1);
                for (let k = 0; k < w * h; k++) {
                    const v = Math.max(0, Math.min(255, 255 * (cinza[k] - lo) / (hi - lo)));
                    dados.data[4 * k] = dados.data[4 * k + 1] = dados.data[4 * k + 2] = v;
                }
                cx.putImageData(dados, 0, 0);
                ok({ w, h, cinza, tela: cv });
            };
            im.onerror = () => { cacheImg.delete(e); falha(new Error(`image ${e} eV`)); };
            im.src = `/iv/${S.sessao}/img/${e}.png`;
        }));
    }
    return cacheImg.get(e);
}

// ------------------------------------------------------------------ calculo (servidor)

async function post(rota, corpo) {
    const resp = await fetch(`/iv/${S.sessao}/${rota}`, {
        method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(corpo),
    });
    if (!resp.ok) {
        let msg = `HTTP ${resp.status}`;
        try { msg = (await resp.json()).erro || msg; } catch (_) { /* resposta sem JSON */ }
        throw new Error(msg);
    }
    return resp;
}

async function recalcula() {
    const n = ++S.pedido;
    S.pos = {}; S.inten = {};
    try {
        const r = await (await post('calcula', {
            ancoras: S.ancoras, metodo: S.metodo, fundo: S.fundo, largura: S.largura, centro: S.centro,
        })).json();
        if (n !== S.pedido) return;  // ja chegou um pedido mais novo
        S.pos = Object.fromEntries(Object.entries(r.pos || {}).map(([e, p]) => [Number(e), p]));
        S.inten = Object.fromEntries(Object.entries(r.inten || {}).map(([e, v]) => [Number(e), v]));
        S.centroAjustado = r.centro_ajustado || null;
        S.mensagem = r.mensagem;
    } catch (erro) {
        if (n !== S.pedido) return;
        S.mensagem = `Error: ${erro.message}`;
    }
    desenha();
}

// ------------------------------------------------------------------ desenho

const roundIdl = v => (v >= 0 ? Math.floor(v + 0.5) : -Math.floor(-v + 0.5));

function preparaCanvas(cv, largCss, altCss) {
    const dpr = window.devicePixelRatio || 1;
    cv.style.height = `${altCss}px`;
    cv.width = Math.round(largCss * dpr);
    cv.height = Math.round(altCss * dpr);
    const cx = cv.getContext('2d');
    cx.setTransform(dpr, 0, 0, dpr, 0, 0);
    cx.imageSmoothingEnabled = false;
    return cx;
}

function energiaAtual() { return S.energias[S.i]; }

function centroAtual() {
    const e = energiaAtual();
    return S.pos[e] || S.ancoras[e] || null;
}

async function desenha() {
    if (!S.sessao) return;
    const e = energiaAtual();
    $('iv-titulo').textContent = `${e} eV`;
    $('iv-slider-valor').textContent = `${e}`;
    desenhaInfo();
    desenhaCurva();
    let im;
    try {
        im = await imagem(e);
    } catch (erro) {
        S.mensagem = `Error loading ${erro.message}.`;
        desenhaInfo();
        return;
    }
    if (energiaAtual() !== e) return;  // o usuario ja passou para outra
    desenhaImagem(im);
    desenhaZoom(im);
}

function desenhaImagem(im) {
    const cv = $('iv-img');
    const larg = cv.clientWidth;
    const s = larg / im.w;
    const cx = preparaCanvas(cv, larg, im.h * s);
    cx.drawImage(im.tela, 0, 0, im.w * s, im.h * s);
    // Pixel k ocupa de k-0.5 a k+0.5 em coordenada de dado, como no imshow.
    const px = v => (v + 0.5) * s;

    const traj = S.energias.filter(k => S.pos[k]);
    if (traj.length) {
        cx.strokeStyle = 'rgba(0, 255, 255, 0.6)';
        cx.lineWidth = 1;
        cx.beginPath();
        traj.forEach((k, j) => (j ? cx.lineTo : cx.moveTo).call(cx, px(S.pos[k][0]), px(S.pos[k][1])));
        cx.stroke();
    }
    const a = S.ancoras[energiaAtual()];
    if (a) {
        cx.strokeStyle = 'yellow';
        cx.lineWidth = 2;
        const x = px(a[0]), y = px(a[1]), r = 5;
        cx.beginPath();
        cx.moveTo(x - r, y - r); cx.lineTo(x + r, y + r);
        cx.moveTo(x - r, y + r); cx.lineTo(x + r, y - r);
        cx.stroke();
    }
    const c = centroAtual();
    if (c) {
        const meio = S.largura / 2;
        cx.strokeStyle = 'lime';
        cx.lineWidth = 1.5;
        cx.strokeRect(px(roundIdl(c[0]) - meio), px(roundIdl(c[1]) - meio), S.largura * s, S.largura * s);
    }
}

function desenhaZoom(im) {
    const cv = $('iv-zoom');
    const larg = cv.clientWidth;
    const cx = preparaCanvas(cv, larg, larg);
    const c = centroAtual();
    if (!c) {
        cx.fillStyle = '#111';
        cx.fillRect(0, 0, larg, larg);
        return;
    }
    const z = 20 + Math.floor(S.largura / 2);
    const n = 2 * z + 1;
    const cc = roundIdl(c[0]), cl = roundIdl(c[1]);
    // Fora da imagem conta como 0, como o np.pad do leed_iv.py.
    const recorte = new Uint8Array(n * n);
    let mn = 255, mx = 0;
    for (let y = 0; y < n; y++) {
        for (let x = 0; x < n; x++) {
            const ix = cc - z + x, iy = cl - z + y;
            const v = ix >= 0 && iy >= 0 && ix < im.w && iy < im.h ? im.cinza[iy * im.w + ix] : 0;
            recorte[y * n + x] = v;
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
    }
    mx = Math.max(mx, mn + 1);
    const tmp = document.createElement('canvas');
    tmp.width = n; tmp.height = n;
    const tx = tmp.getContext('2d');
    const dados = tx.createImageData(n, n);
    for (let k = 0; k < n * n; k++) {
        const v = 255 * (recorte[k] - mn) / (mx - mn);
        dados.data[4 * k] = dados.data[4 * k + 1] = dados.data[4 * k + 2] = v;
        dados.data[4 * k + 3] = 255;
    }
    tx.putImageData(dados, 0, 0);
    const s = larg / n;
    cx.drawImage(tmp, 0, 0, larg, larg);
    const meio = S.largura / 2;
    cx.strokeStyle = 'lime';
    cx.lineWidth = 1.5;
    cx.strokeRect((z - meio + 0.5) * s, (z - meio + 0.5) * s, S.largura * s, S.largura * s);
}

function quebra(texto, n) {
    const linhas = [];
    let atual = '';
    for (const p of texto.split(/\s+/).filter(Boolean)) {
        if (atual.length + p.length + 1 > n) { linhas.push(atual); atual = p; }
        else atual = `${atual} ${p}`.trim();
    }
    if (atual) linhas.push(atual);
    return linhas;
}

function desenhaInfo() {
    const e = energiaAtual();
    const x = S.xml[e] || {};
    const c = centroAtual();
    const l = [`energy      ${e} eV` + (x.Energy ? `  (xml ${x.Energy})` : '')];
    if (x.BeamCurrent) l.push(`current     ${x.BeamCurrent} uA`);
    if (c) l.push(`position    col ${roundIdl(c[0])}  row ${roundIdl(c[1])}`);
    if (e in S.inten) l.push(`intensity   ${S.inten[e]}`);
    const marcadas = Object.keys(S.ancoras).map(Number).sort((a, b) => a - b);
    l.push('', `marked points: ${marcadas.length}`);
    marcadas.forEach(a => l.push(`  ${a} eV: (${S.ancoras[a][0].toFixed(0)}, ${S.ancoras[a][1].toFixed(0)})`));
    l.push('', `window ${S.largura}x${S.largura} px`);
    if (S.centro) l.push(`pattern center (${S.centro[0].toFixed(0)}, ${S.centro[1].toFixed(0)})`);
    l.push('', ...quebra(S.mensagem || '', 38));
    $('iv-info').textContent = l.join('\n');
}

// ------------------------------------------------------------------ curva

function suaviza3(v) {
    if (v.length < 3) return v.slice();
    const p = [v[0], ...v, v.at(-1)];
    return v.map((_, k) => (p[k] + p[k + 1] + p[k + 2]) / 3);
}

const maximo = v => Math.max(...v.filter(Number.isFinite));

function curvaAtual() {
    if (!S.energias.every(e => e in S.inten)) return null;
    let y = S.energias.map(e => S.inten[e]);
    if (S.normalizar && S.temCorrente) y = y.map((v, k) => v / S.xml[S.energias[k]].BeamCurrent);
    if (S.suavizar) y = suaviza3(y);
    const m = maximo(y);
    return y.map(v => v / m);
}

// Mesmo calculo do arquivo salvo (coluna "normalized").
function normalizadaSalva() {
    const inten = S.energias.map(e => S.inten[e]);
    if (S.temCorrente) {
        const pc = inten.map((v, k) => v / S.xml[S.energias[k]].BeamCurrent);
        const m = maximo(pc);
        return pc.map(v => v / m);
    }
    const m = maximo(inten);
    return inten.map(v => v / m);
}

function ticks(min, max, alvo) {
    const passoBruto = (max - min) / alvo;
    const mag = 10 ** Math.floor(Math.log10(passoBruto));
    const passo = [1, 2, 2.5, 5, 10].map(f => f * mag).find(p => p >= passoBruto);
    const out = [];
    for (let v = Math.ceil(min / passo) * passo; v <= max + passo * 1e-9; v += passo) out.push(+v.toFixed(10));
    return out;
}

const geo = { x0: 0, x1: 0, y0: 0, y1: 0, emin: 0, emax: 1, ymin: 0, ymax: 1 };

function desenhaCurva(hoverE = null) {
    const cv = $('iv-curva');
    const larg = cv.clientWidth, alt = 300;
    const cx = preparaCanvas(cv, larg, alt);
    cx.imageSmoothingEnabled = true;
    const es = S.energias;
    if (!es.length) return;

    const series = S.salvos.map((sv, k) => ({
        nome: sv.nome, cor: CORES_SALVOS[k % CORES_SALVOS.length], larg: 1.5, alfa: 0.6,
        y: (() => { const y = S.suavizar ? suaviza3(sv.norm) : sv.norm; const m = maximo(y); return y.map(v => v / m); })(),
    }));
    const atual = curvaAtual();
    if (atual) series.push({ nome: 'current spot', cor: COR_ATUAL, larg: 1.8, alfa: 1, y: atual });

    let ymin = 0, ymax = 1;
    const todos = series.flatMap(s => s.y).filter(Number.isFinite);
    if (todos.length) {
        ymin = Math.min(...todos); ymax = Math.max(...todos);
        const folga = (ymax - ymin || 1) * 0.05;
        ymin -= folga; ymax += folga;
    }
    Object.assign(geo, { x0: 62, x1: larg - 12, y0: 12, y1: alt - 44, emin: es[0], emax: es.at(-1), ymin, ymax });
    const X = e => geo.x0 + (e - geo.emin) / ((geo.emax - geo.emin) || 1) * (geo.x1 - geo.x0);
    const Y = v => geo.y1 - (v - geo.ymin) / ((geo.ymax - geo.ymin) || 1) * (geo.y1 - geo.y0);

    // Grade e eixos recessivos.
    cx.font = '11px "Segoe UI", sans-serif';
    cx.fillStyle = '#666';
    cx.strokeStyle = '#ececec';
    cx.lineWidth = 1;
    cx.textAlign = 'right'; cx.textBaseline = 'middle';
    ticks(ymin, ymax, 5).forEach(v => {
        cx.beginPath(); cx.moveTo(geo.x0, Y(v)); cx.lineTo(geo.x1, Y(v)); cx.stroke();
        cx.fillText(v.toFixed(2).replace(/\.?0+$/, '') || '0', geo.x0 - 6, Y(v));
    });
    cx.textAlign = 'center'; cx.textBaseline = 'top';
    ticks(geo.emin, geo.emax, 8).forEach(v => {
        cx.beginPath(); cx.moveTo(X(v), geo.y1); cx.lineTo(X(v), geo.y1 + 4);
        cx.strokeStyle = '#bbb'; cx.stroke();
        cx.fillText(v, X(v), geo.y1 + 6);
    });
    cx.strokeStyle = '#bbb';
    cx.beginPath(); cx.moveTo(geo.x0, geo.y0); cx.lineTo(geo.x0, geo.y1); cx.lineTo(geo.x1, geo.y1); cx.stroke();
    cx.fillStyle = '#444';
    cx.fillText('energy (eV)', (geo.x0 + geo.x1) / 2, alt - 18);
    cx.save();
    cx.translate(14, (geo.y0 + geo.y1) / 2); cx.rotate(-Math.PI / 2);
    cx.textBaseline = 'middle';
    cx.fillText('intensity' + (S.normalizar && S.temCorrente ? ' / current' : '')
        + (S.suavizar ? ' smoothed' : '') + ' (norm.)', 0, 0);
    cx.restore();

    // Energia atual.
    cx.strokeStyle = COR_ENERGIA; cx.lineWidth = 1;
    cx.beginPath(); cx.moveTo(X(energiaAtual()), geo.y0); cx.lineTo(X(energiaAtual()), geo.y1); cx.stroke();

    series.forEach(s => {
        cx.globalAlpha = s.alfa; cx.strokeStyle = s.cor; cx.lineWidth = s.larg; cx.lineJoin = 'round';
        cx.beginPath();
        let caneta = false;
        s.y.forEach((v, k) => {
            if (!Number.isFinite(v)) { caneta = false; return; }
            caneta ? cx.lineTo(X(es[k]), Y(v)) : cx.moveTo(X(es[k]), Y(v));
            caneta = true;
        });
        cx.stroke();
    });
    cx.globalAlpha = 1;

    // Energias marcadas sobre a curva atual.
    if (atual) {
        cx.strokeStyle = COR_ANCORA; cx.lineWidth = 2;
        Object.keys(S.ancoras).map(Number).forEach(a => {
            const v = atual[es.indexOf(a)];
            if (!Number.isFinite(v)) return;
            const x = X(a), y = Y(v), r = 4;
            cx.beginPath();
            cx.moveTo(x - r, y - r); cx.lineTo(x + r, y + r);
            cx.moveTo(x - r, y + r); cx.lineTo(x + r, y - r);
            cx.stroke();
        });
    }

    // Legenda, so com 2+ curvas (uma so e o "spot atual" e dispensa caixa).
    if (series.length >= 2) {
        cx.font = '11px "Segoe UI", sans-serif';
        const lw = Math.max(...series.map(s => cx.measureText(s.nome).width)) + 34;
        const lx = geo.x1 - lw - 4, ly = geo.y0 + 4, lh = 16;
        cx.fillStyle = 'rgba(255,255,255,0.88)';
        cx.fillRect(lx, ly, lw, lh * series.length + 6);
        cx.textAlign = 'left'; cx.textBaseline = 'middle';
        series.forEach((s, k) => {
            const yy = ly + 11 + k * lh;
            cx.strokeStyle = s.cor; cx.lineWidth = 2;
            cx.beginPath(); cx.moveTo(lx + 6, yy); cx.lineTo(lx + 24, yy); cx.stroke();
            cx.fillStyle = '#333';
            cx.fillText(s.nome, lx + 30, yy);
        });
    }

    // Mira do mouse.
    if (hoverE !== null) {
        cx.strokeStyle = 'rgba(0,0,0,0.35)'; cx.lineWidth = 1; cx.setLineDash([3, 3]);
        cx.beginPath(); cx.moveTo(X(hoverE), geo.y0); cx.lineTo(X(hoverE), geo.y1); cx.stroke();
        cx.setLineDash([]);
        const k = es.indexOf(hoverE);
        series.forEach(s => {
            const v = s.y[k];
            if (!Number.isFinite(v)) return;
            cx.fillStyle = s.cor; cx.strokeStyle = 'white'; cx.lineWidth = 2;
            cx.beginPath(); cx.arc(X(hoverE), Y(v), 4, 0, 2 * Math.PI); cx.fill(); cx.stroke();
        });
        return series.map(s => ({ nome: s.nome, cor: s.cor, v: s.y[k] }));
    }
    return null;
}

function energiaMaisProxima(xCss) {
    const e = geo.emin + (xCss - geo.x0) / ((geo.x1 - geo.x0) || 1) * (geo.emax - geo.emin);
    let melhor = 0;
    S.energias.forEach((k, j) => { if (Math.abs(k - e) < Math.abs(S.energias[melhor] - e)) melhor = j; });
    return melhor;
}

// ------------------------------------------------------------------ eventos

function vaiPara(i) {
    S.i = Math.max(0, Math.min(S.energias.length - 1, i));
    $('iv-slider').value = S.i;
    desenha();
}

async function cliqueImagem(ev) {
    if (!S.sessao) return;
    const cv = $('iv-img');
    const r = cv.getBoundingClientRect();
    const [ny, nx] = S.forma;
    const s = r.width / nx;
    const x = (ev.clientX - r.left) / s - 0.5;
    const y = (ev.clientY - r.top) / s - 0.5;
    if (x < -0.5 || y < -0.5 || x > nx - 0.5 || y > ny - 0.5) return;
    const e = energiaAtual();
    if (ev.button === 0) {
        try {
            S.ancoras[e] = await (await post('marca', { energia: e, x, y })).json();
        } catch (erro) {
            S.mensagem = `Error: ${erro.message}`;
            desenha();
            return;
        }
    } else if (ev.button === 2) {
        delete S.ancoras[e];
    } else {
        return;
    }
    desenha();
    recalcula();
}

async function salva() {
    if (!Object.keys(S.inten).length) {
        S.mensagem = 'Nothing to save: mark the spot first.';
        desenhaInfo();
        return;
    }
    const nome = prompt('Spot name (saves <name>.txt, <name>_pontos.txt and <name>.png in a ZIP):',
        `spot_${S.salvos.length + 1}`);
    if (nome === null) {
        S.mensagem = 'Save cancelled.';
        desenhaInfo();
        return;
    }
    let resp;
    try {
        resp = await post('salva', {
            nome, ancoras: S.ancoras, metodo: S.metodo, fundo: S.fundo, largura: S.largura,
            suavizar: S.suavizar, centro: S.centro,
        });
    } catch (erro) {
        S.mensagem = `Error: ${erro.message}`;
        desenhaInfo();
        return;
    }
    const blob = await resp.blob();
    const arq = (resp.headers.get('Content-Disposition') || '').match(/filename="?([^";]+)"?/);
    const nomeZip = arq ? arq[1] : `${nome}.zip`;
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = nomeZip;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 10000);

    if (S.metodo === 'fisico' && Object.keys(S.ancoras).length >= 2 && S.centroAjustado) {
        S.centro = S.centroAjustado;
    }
    S.salvos.push({ nome: nomeZip.replace(/\.zip$/, ''), norm: normalizadaSalva() });
    S.mensagem = `Saved ${nomeZip}. 'New spot' to measure another one.`;
    desenha();
}

function liga() {
    $('iv-pasta-img').addEventListener('change', atualizaContagem);
    $('iv-pasta-xml').addEventListener('change', atualizaContagem);
    $('iv-carrega').addEventListener('click', carrega);

    const cvImg = $('iv-img');
    cvImg.addEventListener('mousedown', cliqueImagem);
    cvImg.addEventListener('contextmenu', ev => ev.preventDefault());

    const roda = ev => {
        if (!S.sessao) return;
        ev.preventDefault();
        vaiPara(S.i + (ev.deltaY < 0 ? 1 : -1));
    };
    [cvImg, $('iv-zoom'), $('iv-curva')].forEach(el => el.addEventListener('wheel', roda, { passive: false }));

    const cvCurva = $('iv-curva');
    cvCurva.addEventListener('click', ev => {
        if (!S.sessao) return;
        vaiPara(energiaMaisProxima(ev.clientX - cvCurva.getBoundingClientRect().left));
    });
    cvCurva.addEventListener('mousemove', ev => {
        if (!S.sessao) return;
        const r = cvCurva.getBoundingClientRect();
        const xc = ev.clientX - r.left;
        const tip = $('iv-tooltip');
        if (xc < geo.x0 || xc > geo.x1) { tip.hidden = true; desenhaCurva(); return; }
        const e = S.energias[energiaMaisProxima(xc)];
        const vals = desenhaCurva(e) || [];
        const linhas = vals.filter(v => Number.isFinite(v.v))
            .map(v => `<div><span style="display:inline-block;width:10px;height:2px;background:${v.cor};vertical-align:middle;margin-right:6px"></span>${v.nome}: ${v.v.toFixed(3)}</div>`);
        tip.innerHTML = `<div><b>${e} eV</b></div>${linhas.join('')}`;
        tip.hidden = false;
        const esquerda = xc > r.width / 2 ? xc - tip.offsetWidth - 12 : xc + 12;
        tip.style.left = `${esquerda}px`;
        tip.style.top = `${geo.y0 + 8}px`;
    });
    cvCurva.addEventListener('mouseleave', () => { $('iv-tooltip').hidden = true; desenhaCurva(); });

    $('iv-slider').addEventListener('input', ev => vaiPara(Number(ev.target.value)));

    document.addEventListener('keydown', ev => {
        if (!S.sessao || ev.target.matches('input[type=text], input[type=number], input[type=range], textarea')) return;
        if (ev.key === 'ArrowRight' || ev.key === 'ArrowUp') { ev.preventDefault(); vaiPara(S.i + 1); }
        else if (ev.key === 'ArrowLeft' || ev.key === 'ArrowDown') { ev.preventDefault(); vaiPara(S.i - 1); }
        else if (ev.key === 'Delete' || ev.key === 'Backspace') {
            delete S.ancoras[energiaAtual()];
            recalcula();
        }
    });

    document.querySelectorAll('input[name=iv-fundo]').forEach(r => r.addEventListener('change', () => {
        S.fundo = r.value; recalcula();
    }));
    document.querySelectorAll('input[name=iv-metodo]').forEach(r => r.addEventListener('change', () => {
        S.metodo = r.value; recalcula();
    }));
    $('iv-ck-corrente').addEventListener('change', ev => { S.normalizar = ev.target.checked; desenha(); });
    $('iv-ck-suave').addEventListener('change', ev => { S.suavizar = ev.target.checked; desenha(); });

    $('iv-limpa').addEventListener('click', () => { S.ancoras = {}; recalcula(); });
    $('iv-novo').addEventListener('click', async () => {
        S.ancoras = {};
        await recalcula();
        if (S.centro) {
            S.mensagem = 'New spot: with the pattern center already known, 1 click is enough.';
            desenhaInfo();
        }
    });
    $('iv-salva').addEventListener('click', salva);

    let espera;
    window.addEventListener('resize', () => { clearTimeout(espera); espera = setTimeout(desenha, 100); });
}

liga();
