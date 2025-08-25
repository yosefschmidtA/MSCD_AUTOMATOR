# tradutor_input.py
# Este script lê um arquivo de entrada e gera um cabeçalho formatado para o MSCD.

atomic_weights = {
    # Período 1
    1: 1.008,    # Hidrogênio (H)
    2: 4.0026,   # Hélio (He)
    
    # Período 2
    3: 6.94,     # Lítio (Li)
    4: 9.0122,   # Berílio (Be)
    5: 10.81,    # Boro (B)
    6: 12.011,   # Carbono (C)
    7: 14.007,   # Nitrogênio (N)
    8: 15.999,   # Oxigênio (O)
    9: 18.998,   # Flúor (F)
    10: 20.180,  # Neônio (Ne)
    
    # Período 3
    11: 22.990,  # Sódio (Na)
    12: 24.305,  # Magnésio (Mg)
    13: 26.982,  # Alumínio (Al)
    14: 28.085,  # Silício (Si)
    15: 30.974,  # Fósforo (P)
    16: 32.06,   # Enxofre (S)
    17: 35.45,   # Cloro (Cl)
    18: 39.948,  # Argônio (Ar)
    
    # Período 4
    19: 39.098,  # Potássio (K)
    20: 40.078,  # Cálcio (Ca)
    21: 44.956,  # Escândio (Sc)
    22: 47.867,  # Titânio (Ti)
    23: 50.942,  # Vanádio (V)
    24: 51.996,  # Cromo (Cr)
    25: 54.938,  # Manganês (Mn)
    26: 55.845,  # Ferro (Fe)
    27: 58.933,  # Cobalto (Co)
    28: 58.693,  # Níquel (Ni)
    29: 63.546,  # Cobre (Cu)
    30: 65.38,   # Zinco (Zn)
    31: 69.723,  # Gálio (Ga)
    32: 72.630,  # Germânio (Ge)
    33: 74.922,  # Arsênio (As)
    34: 78.971,  # Selênio (Se)
    35: 79.904,  # Bromo (Br)
    36: 83.798,  # Criptônio (Kr)
    
    # Outros elementos comuns
    47: 107.868, # Prata (Ag)
    74: 183.84,  # Tungstênio (W)
    78: 195.08,  # Platina (Pt)
    79: 196.97,  # Ouro (Au)
    80: 200.59,  # Mercúrio (Hg)
}

def gerar_cabecalho(input_file, system_name="Ga2O3(100)", user="Yosef (UFJ-PPGCAS)"):
    """
    Gera um cabeçalho formatado para um arquivo de entrada de simulação.
    """
    elementos = []
    emissor = None
    exp_file = None
    out_file = None

    try:
        with open(input_file, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        return f"ERRO: O ARQUIVO '{input_file}' NAO FOI ENCONTRADO."

    # === Processar as 4 primeiras linhas ===
    
    # A. Lê a primeira linha e extrai os arquivos exp e out
    first_line_parts = lines[0].split()
    if len(first_line_parts) >= 6:
        exp_file = first_line_parts[4]
        out_file = first_line_parts[5]
    
    # B. Procura por todos os elementos e o emissor
    for i in range(4):
        parts = lines[i].split()
        if not parts:
            continue
            
        # Adiciona o elemento à lista de elementos
        if len(parts) >= 2 and parts[1].isdigit():
            z = int(parts[1])
            if z not in elementos:
                elementos.append(z)
        
        # Lógica corrigida para encontrar o emissor
        if len(parts) >= 7 and parts[-1].strip() == "1":
            emissor = int(parts[1])
        elif i > 0 and parts[-1].strip() == "1":
            emissor = int(parts[1])

    # === Construção do cabeçalho ===
    header = []
    header.append("    741   10    191     datakind begining-row linenumbers\n")
    header.append("----------------------------------------------------------------\n")
    header.append("MSCD Version 1.00 Yufeng Chen and Michel A Van Hove\n")
    header.append("Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720\n")
    header.append("Copyright (c) Van Hove Group 1997. All rights reserved\n")
    header.append("----------------------------------------------------------------\n\n")
    header.append(f"        {system_name}              input file\n\n")
    header.append(f'un      "{user}"    user name\n')
    header.append(f"sn      {system_name}              system name\n")

    for i, z in enumerate(elementos, start=1):
        header.append(f"ps{i:02d}    ps{i}.1.txt               input phase shift data file\n")

    header.append("\n")

    if emissor:
        idx = elementos.index(emissor) + 1
        header.append(f"rm{idx:02d}    rm{idx}.txt               input radial matrix data file\n")

    if exp_file and out_file:
        header.append(f"ex      {exp_file:<20} experimental data\n")
        header.append(f"pe      {out_file:<20} output photoemission data file\n")

    header.append("\n")

    # === Agora copiar SOMENTE até a linha "fit try ..." ===
    for line in lines[4:]:
        header.append(line)
        if "fit try for vinner" in line:
            break

    return "".join(header)

if __name__ == "__main__":
    cabecalho = gerar_cabecalho("input_cluster.txt")
    
    if cabecalho.startswith("ERRO"):
        print(cabecalho)
    else:
        print(cabecalho)
        with open("output_header.txt", "w") as f:
            f.write(cabecalho)
