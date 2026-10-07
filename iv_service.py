"""Leitor de curvas IV de LEED no site: o leed_iv.py do iv_curves, com a janela no navegador.

O navegador manda as duas pastas uma vez; daí em diante cada pedido traz so o estado
da tela (pontos marcados, metodo, fundo) e o calculo roda aqui, no leed_iv_core.py,
que reproduz o coleta_no_step do IDL.
"""
import io
import re
import shutil
import tempfile
import time
import uuid
from pathlib import Path

import numpy as np
from PIL import Image

import leed_iv_core as core

# Fora da pasta do projeto: no dev o codigo e montado do host e as imagens enviadas
# iriam parar no git. Somem quando o conteiner e recriado, e e isso que se quer.
RAIZ = Path(tempfile.gettempdir()) / "iv_sessoes"
VALIDADE_S = 24 * 3600

NOME_IMAGEM = re.compile(r"^(\d+)\.(jpe?g|tiff?)$", re.IGNORECASE)
NOME_XML = re.compile(r"^(\d+)\.xml$", re.IGNORECASE)
NOME_SESSAO = re.compile(r"^[0-9a-f]{32}$")


class ErroIV(Exception):
    pass


def _limpa_antigas():
    if not RAIZ.exists():
        return
    agora = time.time()
    for d in RAIZ.iterdir():
        if d.is_dir() and agora - d.stat().st_mtime > VALIDADE_S:
            shutil.rmtree(d, ignore_errors=True)


def _pasta(sessao):
    if not NOME_SESSAO.fullmatch(sessao or ""):
        raise ErroIV("Invalid session.")
    d = RAIZ / sessao
    if not d.is_dir():
        raise ErroIV("Session expired or not found. Upload the folders again.")
    # A validade conta do ultimo uso, nao do envio: quem ainda esta medindo nao
    # perde a sessao para a limpeza disparada pelo Load de outra pessoa.
    d.touch()
    return d


def cria_sessao(arquivos_img, arquivos_xml):
    """Grava as imagens <energia>.jpg/.tif e os <energia>.xml; o resto e ignorado.

    Os XML podem vir na pasta propria ou junto das imagens, como no leed_iv.py.
    """
    _limpa_antigas()
    sessao = uuid.uuid4().hex
    d = RAIZ / sessao
    (d / "img").mkdir(parents=True)
    (d / "xml").mkdir()
    for f in list(arquivos_img) + list(arquivos_xml):
        # Com a pasta inteira o navegador manda "pasta/100.jpg"; so o nome interessa.
        nome = Path(f.filename or "").name
        if NOME_IMAGEM.fullmatch(nome):
            f.save(d / "img" / nome)
        elif NOME_XML.fullmatch(nome):
            f.save(d / "xml" / nome)
    if not core.lista_imagens(d / "img"):
        shutil.rmtree(d, ignore_errors=True)
        raise ErroIV("No images named <energy>.jpg / .tif (e.g. 100.jpg) in the images folder.")
    return sessao


def _imagens(d):
    return core.lista_imagens(d / "img")


def _xml(d, energias):
    # A energia vem do nome do arquivo; do XML so a corrente entra no calculo.
    return {e: core.le_xml(d / "xml" / f"{e}.xml") for e in energias}


def info(sessao):
    d = _pasta(sessao)
    imagens = _imagens(d)
    energias = list(imagens)
    xml = _xml(d, energias)
    ny, nx = core.carrega_cinza(imagens[energias[0]]).shape
    return {
        "sessao": sessao,
        "energias": energias,
        "largura": int(round(15 * nx / 640)) | 1,  # 15 px do IDL para 640 px
        "forma": [ny, nx],
        "xml": {e: xml[e] for e in energias},
        "tem_corrente": all(v["BeamCurrent"] for v in xml.values()),
        "sem_xml": [e for e in energias if xml[e]["BeamCurrent"] is None],
    }


def png_cinza(sessao, energia):
    """A imagem em cinza, sem contraste: o navegador estica so a tela, como o tvscl."""
    d = _pasta(sessao)
    imagens = _imagens(d)
    if energia not in imagens:
        raise ErroIV(f"No image for {energia} eV.")
    buf = io.BytesIO()
    Image.fromarray(np.asarray(core.carrega_cinza(imagens[energia]))).save(buf, "PNG")
    buf.seek(0)
    return buf


def _ancoras(bruto):
    return {int(e): (float(c), float(l)) for e, (c, l) in (bruto or {}).items()}


def _faixa(imagens, emin, emax):
    """So as energias de emin a emax (cada spot tem a sua faixa), como o --emin/--emax do leed_iv.py."""
    faixa = {e: p for e, p in imagens.items()
             if (emin is None or e >= emin) and (emax is None or e <= emax)}
    if not faixa:
        raise ErroIV(f"No images between {emin} and {emax} eV.")
    return faixa


def marca(sessao, energia, x, y):
    """O clique a mao erra 1-2 px; puxa para o pico mais proximo se houver um claro."""
    d = _pasta(sessao)
    img = core.carrega_cinza(_imagens(d)[energia])
    c, l, _ = core.pico_local(img, x, y, raio=3)
    return [float(core.round_idl(c)), float(core.round_idl(l))]


def calcula(sessao, ancoras, metodo, fundo, largura, centro=None, emin=None, emax=None):
    """Trajetoria e intensidade; mesmas regras e mensagens do Coleta.recalcula."""
    d = _pasta(sessao)
    # A janela agora vem da tela (botoes - / +); impar para ter um pixel central.
    if not (3 <= largura <= 99 and largura % 2):
        raise ErroIV("The window must be an odd number of pixels from 3 to 99.")
    imagens = _faixa(_imagens(d), emin, emax)
    # Marcacao fora da faixa nao entra no ajuste (o navegador ja as descarta).
    ancoras = {e: a for e, a in _ancoras(ancoras).items() if e in imagens}
    centro = tuple(centro) if centro else None
    if not ancoras:
        return {"mensagem": "Click on the center of the spot you want to measure."}
    if metodo == "polinomio" and len(ancoras) < 4:
        return {"mensagem": f"Polynomial (as in IDL) needs 4 points or more; there are {len(ancoras)}."}
    if metodo == "fisico" and len(ancoras) == 1 and centro is None:
        return {"mensagem": "Mark the same spot at another energy, far from this one."}
    try:
        pos = core.trajetoria(imagens, ancoras, metodo, largura, centro=centro)
    except ValueError as erro:
        return {"mensagem": str(erro)}
    inten = core.curva(imagens, pos, fundo, largura)
    resposta = {
        "pos": {e: list(p) for e, p in pos.items()},
        "inten": inten,
        "mensagem": (f"{len(ancoras)} point(s) marked. Check the green square along the "
                     "energies; if it leaves the spot, click on the spot there."),
    }
    if metodo == "fisico" and len(ancoras) >= 2:
        # Depois de gravar, o navegador guarda isto e o proximo spot sai com 1 clique.
        cc, _, cl, _ = core.ajusta_fisico(ancoras)
        resposta["centro_ajustado"] = [float(cc), float(cl)]
    return resposta


def salva(sessao, nome, ancoras, metodo, fundo, largura, suavizar, centro=None, emin=None, emax=None):
    """So o <nome>.txt da curva, no formato do leed_iv.py (os pontos marcados vao no cabecalho)."""
    d = _pasta(sessao)
    imagens = _faixa(_imagens(d), emin, emax)
    energias = list(imagens)
    # Corrente e normalizacao so dentro da faixa, como no leed_iv.py com --emin/--emax.
    xml = _xml(d, energias)
    tem_corrente = all(v["BeamCurrent"] for v in xml.values())
    ancoras = {e: a for e, a in _ancoras(ancoras).items() if e in imagens}
    nome = re.sub(r"[^\w.-]", "_", nome or "spot_1").removesuffix(".txt") or "spot_1"

    r = calcula(sessao, ancoras, metodo, fundo, largura, centro, emin, emax)
    if "pos" not in r:
        raise ErroIV("Nothing to save: mark the spot first.")
    pos, inten_d = r["pos"], r["inten"]

    es = energias
    inten = np.array([inten_d[e] for e in es], dtype=float)
    suave = core.suaviza3(inten)
    corr = np.array([xml[e]["BeamCurrent"] or np.nan for e in es], dtype=float)
    por_corr = inten / corr
    norm = por_corr / np.nanmax(por_corr) if tem_corrente else inten / inten.max()
    norm_suave = core.suaviza3(norm)

    txt = io.StringIO()
    txt.write(f"# LEED IV curve - {nome}\n")
    txt.write(f"# position method: {metodo} | background: {fundo} | window: {largura} px\n")
    txt.write(f"# energy range: {energias[0]} to {energias[-1]} eV\n")
    txt.write("# marked points (energy: column, row): "
              + "; ".join(f"{a}: {ancoras[a][0]:.0f}, {ancoras[a][1]:.0f}" for a in sorted(ancoras)) + "\n")
    txt.write("# intensity and smoothed: same calculation as coleta_no_step (IDL), without losing the last energy\n")
    txt.write("# normalized = intensity / current, divided by the maximum; normalized_smooth = 3-point average\n")
    txt.write(f"# smoothing on screen when saved: {'yes' if suavizar else 'no'}\n")
    txt.write("# energy energy_xml column row intensity smoothed current_uA intensity_per_uA normalized normalized_smooth\n")
    for k, e in enumerate(es):
        ex = xml[e]["Energy"]
        txt.write(f"{e:5d} {ex if ex is not None else float('nan'):7.1f} {pos[e][0]:6.0f} {pos[e][1]:6.0f} "
                  f"{inten[k]:10.0f} {suave[k]:12.4f} {corr[k]:6.2f} {por_corr[k]:12.2f} "
                  f"{norm[k]:10.5f} {norm_suave[k]:10.5f}\n")

    buf = io.BytesIO(txt.getvalue().encode("utf-8"))
    return nome, buf
