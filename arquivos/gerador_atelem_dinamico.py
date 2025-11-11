import os
import subprocess
import time

# --- Dicionário de Símbolos Atômicos (completo) ---
ATOMIC_SYMBOLS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt',
    110: 'Ds', 111: 'Rg', 112: 'Cn'
}

# --- LÓGICA DE GERAÇÃO DA CONFIGURAÇÃO ELETRÔNICA ---
AUFBAU_ORDER = [(1,0),(2,0),(2,1),(3,0),(3,1),(4,0),(3,2),(4,1),(5,0),(4,2),(5,1),(6,0),(4,3),(5,2),(6,1),(7,0),(5,3),(6,2),(7,1)]
EXCEPTIONS = {24:('Cr',(4,0),(3,2)), 29:('Cu',(4,0),(3,2)), 41:('Nb',(5,0),(4,2)), 42:('Mo',(5,0),(4,2)), 44:('Ru',(5,0),(4,2)), 45:('Rh',(5,0),(4,2)), 46:('Pd',(5,0),(4,2)), 47:('Ag',(5,0),(4,2)), 78:('Pt',(6,0),(5,2)), 79:('Au',(6,0),(5,2))}
def get_electronic_configuration(z):
    electrons_to_place=z; config={};
    for n,l in AUFBAU_ORDER:
        if electrons_to_place==0: break
        capacity=2*(2*l+1); electrons=min(electrons_to_place,capacity); config[(n,l)]=electrons; electrons_to_place-=electrons
    if z in EXCEPTIONS:
        symbol,s_orbital,d_orbital=EXCEPTIONS[z]
        if s_orbital in config and d_orbital in config and config[s_orbital]>0:
            electrons_to_move=1 if z!=46 else 2
            if config[s_orbital]>=electrons_to_move: config[s_orbital]-=electrons_to_move; config[d_orbital]+=electrons_to_move
    final_config=[]
    for(n,l),electrons in sorted(config.items()):
        if electrons==0: continue
        if l==0: final_config.append((n,l,0.5,electrons))
        else:
            j_minus,j_plus=l-0.5,l+0.5; capacity_j_minus=2*j_minus+1
            electrons_in_j_minus=min(electrons,capacity_j_minus); final_config.append((n,l,j_minus,electrons_in_j_minus))
            electrons_left=electrons-electrons_in_j_minus
            if electrons_left>0: final_config.append((n,l,j_plus,electrons_left))
    return final_config

def create_and_run(atomic_number):
    symbol = ATOMIC_SYMBOLS.get(atomic_number)
    if not symbol: return

    print(f"\n--- Processando Elemento: {symbol} (Z={atomic_number}) ---")

    config = get_electronic_configuration(atomic_number)
    atorb_filename = f"Atorb.{symbol}"
    output_filename = f"atelem.{symbol}"
    
    # Parâmetros de convergência iniciais
    mixing, etol, xnum = 0.3, 0.00001, 300
    
    for attempt in range(1, 4):
        print(f"  -> Tentativa {attempt} com mixing={mixing:.1f}, etol={etol:.1e}, xnum={xnum}")

        num_levels = len(config)
        content = ["i", f"{atomic_number} 1000", "d", "1", "x", "0.0d0", "a", f"0 {num_levels} {mixing} {etol} {xnum}"]
        for n, l, j, occ in config: content.append(f"{n} {l} {l} {-j} 1  {occ:.8f}")
        content.extend(["w", output_filename, "q"])
        
        with open(atorb_filename, 'w') as f: f.write("\n".join(content) + "\n")
        print(f"  -> Arquivo '{atorb_filename}' criado.")

        executable = './phsh0'
        if not os.path.exists(executable):
            print(f"ERRO: Executável '{executable}' não encontrado."); return
            
        start_time = time.strftime('%H:%M:%S')
        print(f"  -> [{start_time}] Executando (Timeout: 60s)...")
        try:
            with open(atorb_filename, 'r') as stdin_file:
                subprocess.run([executable], stdin=stdin_file, text=True, check=True, timeout=60)
            
            end_time = time.strftime('%H:%M:%S')
            print(f"  -> [{end_time}] SUCESSO! Arquivo '{output_filename}' gerado.")
            return

        except subprocess.TimeoutExpired:
            end_time = time.strftime('%H:%M:%S')
            print(f"  -> [{end_time}] AVISO: Cálculo excedeu 60 segundos.")
            print("     Modificando parâmetros para a próxima tentativa.")
            # SUA NOVA LÓGICA DE CONVERGÊNCIA
            etol *= 100
            mixing += 0.2
            xnum += 50

        except subprocess.CalledProcessError as e:
            print(f"  ERRO durante a execução do phsh0 para {symbol}:")
            if e.stderr: print(f"--- Saída de Erro ---\n{e.stderr}\n---------------------")
            return

    print(f"  FALHA: O cálculo para {symbol} não convergiu após 3 tentativas.")

# --- BLOCO PRINCIPAL ---
if __name__ == "__main__":
    atomic_numbers_to_run = range(25, 90)
    print("="*50)
    print("INICIANDO GERADOR DINÂMICO DE ARQUIVOS ATELEM")
    print("="*50)
    for z in atomic_numbers_to_run:
        create_and_run(z)
        time.sleep(1)
    print("\nProcesso concluído.")
