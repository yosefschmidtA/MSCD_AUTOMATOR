# input_translator.py
# This script reads an input file and generates a formatted header for MSCD.

atomic_weights = {
    # Period 1
    1: 1.008,    # Hydrogen (H)
    2: 4.0026,   # Helium (He)
    
    # Period 2
    3: 6.94,     # Lithium (Li)
    4: 9.0122,   # Beryllium (Be)
    5: 10.81,    # Boron (B)
    6: 12.011,   # Carbon (C)
    7: 14.007,   # Nitrogen (N)
    8: 15.999,   # Oxygen (O)
    9: 18.998,   # Fluorine (F)
    10: 20.180,  # Neon (Ne)
    
    # Period 3
    11: 22.990,  # Sodium (Na)
    12: 24.305,  # Magnesium (Mg)
    13: 26.982,  # Aluminum (Al)
    14: 28.085,  # Silicon (Si)
    15: 30.974,  # Phosphorus (P)
    16: 32.06,   # Sulfur (S)
    17: 35.45,   # Chlorine (Cl)
    18: 39.948,  # Argon (Ar)
    
    # Period 4
    19: 39.098,  # Potassium (K)
    20: 40.078,  # Calcium (Ca)
    21: 44.956,  # Scandium (Sc)
    22: 47.867,  # Titanium (Ti)
    23: 50.942,  # Vanadium (V)
    24: 51.996,  # Chromium (Cr)
    25: 54.938,  # Manganese (Mn)
    26: 55.845,  # Iron (Fe)
    27: 58.933,  # Cobalt (Co)
    28: 58.693,  # Nickel (Ni)
    29: 63.546,  # Copper (Cu)
    30: 65.38,   # Zinc (Zn)
    31: 69.723,  # Gallium (Ga)
    32: 72.630,  # Germanium (Ge)
    33: 74.922,  # Arsenic (As)
    34: 78.971,  # Selenium (Se)
    35: 79.904,  # Bromine (Br)
    36: 83.798,  # Krypton (Kr)
    
    # Other common elements
    47: 107.868, # Silver (Ag)
    74: 183.84,  # Tungsten (W)
    78: 195.08,  # Platinum (Pt)
    79: 196.97,  # Gold (Au)
    80: 200.59,  # Mercury (Hg)
}

def gerar_cabecalho(input_file, system_name="Ga2O3(100)", user="Yosef (UFJ-PPGCAS)"):
    """
    Generates a formatted header for a simulation input file.
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

    # === Process the first 4 lines ===
    
    # A. Reads the first line and extracts the exp and out files
    first_line_parts = lines[0].split()
    if len(first_line_parts) >= 6:
        exp_file = first_line_parts[4]
        out_file = first_line_parts[5]
    
    # B. Searches for all elements and the emitter
    for i in range(4):
        parts = lines[i].split()
        if not parts:
            continue
            
        # Adds the element to the element list
        if len(parts) >= 2 and parts[1].isdigit():
            z = int(parts[1])
            if z not in elementos:
                elementos.append(z)
        
        # Corrected logic to find the emitter
        if len(parts) >= 7 and parts[-1].strip() == "1":
            emissor = int(parts[1])
        elif i > 0 and parts[-1].strip() == "1":
            emissor = int(parts[1])

    # === Header construction ===
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

    # === Now copy ONLY up to the "fit try ..." line ===
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
