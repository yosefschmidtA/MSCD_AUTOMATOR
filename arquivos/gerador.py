
import numpy as np
import math
import re
import os
import subprocess

# --- Constantes de Conversão ---
A_TO_BOHR = 1.88973

# Dicionário para converter NÚMERO ATÔMICO para SÍMBOLO
atomic_symbols = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 47: 'Ag',
    74: 'W'
}

# --- NOVA FUNÇÃO ---
def get_lattice_constant_from_header(filename='input_cluster.txt'):
    """Lê o cabeçalho para encontrar o parâmetro de rede (lattice)."""
    print(f"\nProcurando pelo parâmetro de rede (lattice) em: {filename}")
    try:
        with open(filename, 'r') as f:
            for line in f:
                # A linha de interesse contém 'lattice(angs)'
                if 'lattice(angs)' in line:
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        lattice_val = float(parts[2])
                        print(f"Parâmetro de rede encontrado: {lattice_val} Angstroms. Este valor será usado como multiplicador.")
                        return lattice_val
                    else:
                        print("AVISO: Linha 'lattice(angs)' encontrada, mas malformada. Usando padrão 1.0.")
                        return 1.0
                # Para de procurar quando os blocos de camada começam
                if '///' in line:
                    break
        print("AVISO: Parâmetro de rede 'lattice(angs)' não encontrado no cabeçalho. Usando padrão 1.0 (unidades em Angstroms).")
        return 1.0
    except Exception as e:
        print(f"ERRO ao ler o parâmetro de rede do cabeçalho: {e}. Usando padrão 1.0.")
        return 1.0

# --- FUNÇÃO MODIFICADA ---
def parse_multi_layer_file(filename='input_cluster.txt', lattice_constant=1.0):
    """Lê um arquivo com múltiplas layers, aplicando o fator de rede aos comprimentos."""
    print(f"Lendo o arquivo de configuração do cluster: {filename}")
    try:
        with open(filename, 'r') as f:
            all_lines = f.readlines()

        first_layer_start_index = -1
        for i, line in enumerate(all_lines):
            if 'unita' in line:
                for j in range(i - 1, -1, -1):
                    parts = all_lines[j].strip().split()
                    if len(parts) >= 2:
                        try:
                            int(parts[0])
                            int(parts[1])
                            first_layer_start_index = j
                            break
                        except (ValueError, IndexError):
                            continue
                if first_layer_start_index != -1:
                    break
        
        if first_layer_start_index == -1:
            print("ERRO: Não foi possível encontrar o início dos dados da camada (marcador 'unita').")
            return None

        content = "".join(all_lines[first_layer_start_index:])
        layer_blocks = re.split(r'/{3,}', content)
        all_layers_params = []

        print(f"Aplicando fator de rede de {lattice_constant} aos comprimentos.")
        for block in layer_blocks:
            if not block.strip(): continue
            params = {}
            lines = [line for line in block.strip().split('\n') if line.strip()]
            parts = lines[0].split()
            params['atomic_number'] = int(parts[0])
            params['is_emitter'] = (int(parts[1]) == 1)
            for i, line in enumerate(lines):
                if 'unita' in line:
                    val_len = float(lines[i + 1].split()[0])
                    val_ang = float(lines[i + 1].split()[1])
                    params.update({'unita_len': val_len * lattice_constant, 'unita_ang_deg': val_ang})
                elif 'unitb' in line:
                    val_len = float(lines[i + 1].split()[0])
                    val_ang = float(lines[i + 1].split()[1])
                    params.update({'unitb_len': val_len * lattice_constant, 'unitb_ang_deg': val_ang})
                elif 'origin' in line:
                    val_len = float(lines[i + 1].split()[0])
                    val_ang = float(lines[i + 1].split()[1])
                    params.update({'origin_len': val_len * lattice_constant, 'origin_ang_deg': val_ang})
                elif 'interlayer' in line:
                    val_spacing = float(lines[i + 1].split()[0])
                    params.update({'interlayer_spacing': val_spacing * lattice_constant})
            all_layers_params.append(params)
        print(f"Encontradas {len(all_layers_params)} layers para processar.")
        return all_layers_params
    except Exception as e:
        print(f"ERRO ao ler o arquivo de input: {e}")
        return None

def get_element_order_from_header(filename='input_cluster.txt'):
    """Lê o cabeçalho do arquivo para obter a ordem dos elementos para o mufftin.d."""
    print(f"\nLendo a ordem dos elementos do cabeçalho de: {filename}")
    ordered_atomic_numbers = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split()
                try:
                    int(parts[0])
                    atomic_num = int(parts[1])
                    float(parts[3])
                    if atomic_num not in ordered_atomic_numbers:
                        ordered_atomic_numbers.append(atomic_num)
                except (ValueError, IndexError):
                    break
        print(f"Ordem de elementos para mufftin.d definida pelo cabeçalho: {ordered_atomic_numbers}")
        return ordered_atomic_numbers
    except Exception as e:
        print(f"ERRO ao ler o cabeçalho do arquivo de input: {e}")
        return []


def assemble_mufftin_d(ordered_atomic_numbers, output_filename='mufftin.d'):
    """Monta o arquivo mufftin.d a partir dos arquivos atelem.*."""
    print(f"\nMontando o arquivo {output_filename} a partir dos arquivos atelem.*...")
    try:
        if not ordered_atomic_numbers:
            print("Nenhuma ordem de elemento fornecida para montar o arquivo.")
            return

        symbols = [atomic_symbols.get(z, None) for z in ordered_atomic_numbers]
        symbols = [s for s in symbols if s is not None]
        print(f"Ordem dos elementos para montagem: {', '.join(symbols)}")

        atelem_files = []
        for symbol in symbols:
            filename = f"atelem.{symbol}"
            if os.path.exists(filename):
                atelem_files.append(filename)
            else:
                print(f"AVISO: Arquivo de potencial '{filename}' não encontrado. Ele será pulado.")

        if not atelem_files:
            print("Nenhum arquivo atelem.* correspondente foi encontrado. O arquivo mufftin.d não será criado.")
            return

        first_file = atelem_files.pop(0)
        command_cp = ['cp', first_file, output_filename]
        print(f"-> Executando: {' '.join(command_cp)}")
        subprocess.run(command_cp, check=True)

        for file_to_append in atelem_files:
            command_cat = f"cat {file_to_append} >> {output_filename}"
            print(f"-> Executando: {command_cat}")
            subprocess.run(command_cat, shell=True, check=True)
        
        print(f"Arquivo '{output_filename}' montado com sucesso.")

    except Exception as e:
        print(f"ERRO ao montar o arquivo {output_filename}: {e}")


def generate_mscd_phaseshift_input(all_params, output_filename='mscd_phaseshift.inp'):
    """Gera o arquivo de input para o MSCD com vetores de rede totalmente adaptativos."""
    print(f"\nGerando arquivo de input para o MSCD: {output_filename}")

    first_layer = all_params[0]
    lattice_constant_A = min(l['unita_len'] for l in all_params if l['unita_len'] > 0)
    lattice_constant_au = lattice_constant_A * A_TO_BOHR

    unita_ang_rad = np.radians(first_layer['unita_ang_deg'])
    unitb_ang_rad = np.radians(first_layer['unitb_ang_deg'])

    unita_cartesian = np.array([
        first_layer['unita_len'] * np.cos(unita_ang_rad),
        first_layer['unita_len'] * np.sin(unita_ang_rad)
    ])
    unitb_cartesian = np.array([
        first_layer['unitb_len'] * np.cos(unitb_ang_rad),
        first_layer['unitb_len'] * np.sin(unitb_ang_rad)
    ])

    vec_a = np.array([unita_cartesian[0] / lattice_constant_A, unita_cartesian[1] / lattice_constant_A, 0.0])
    vec_b = np.array([unitb_cartesian[0] / lattice_constant_A, unitb_cartesian[1] / lattice_constant_A, 0.0])
    vec_c = np.array([0.0, 0.0, 10.0])

    print("Parâmetros de rede calculados (agora baseados em Angstroms):")
    print(f"  - Parâmetro de rede base (a_lat): {lattice_constant_A:.4f} Å")
    print(f"  - Vetor A (normalizado): ({vec_a[0]:.4f}, {vec_a[1]:.4f}, {vec_a[2]:.4f})")
    print(f"  - Vetor B (normalizado): ({vec_b[0]:.4f}, {vec_b[1]:.4f}, {vec_b[2]:.4f})")

    atom_types = {}
    accumulated_z = 0.0
    for params in all_params:
        z_number = params['atomic_number']
        if z_number not in atom_types:
            atom_types[z_number] = []
        origin_ang_rad = np.radians(params['origin_ang_deg'])
        origin_x = params['origin_len'] * np.cos(origin_ang_rad)
        origin_y = params['origin_len'] * np.sin(origin_ang_rad)
        accumulated_z -= params.get('interlayer_spacing', 0.0) # Usar .get para segurança
        
        basis_atom = np.array([
            origin_x / lattice_constant_A,
            origin_y / lattice_constant_A,
            accumulated_z / lattice_constant_A
        ])
        atom_types[z_number].append(basis_atom)

    with open(output_filename, 'w') as f:
        f.write("Input gerado automaticamente para MSCD (\n")
        f.write(f"  {lattice_constant_au:<8.4f}              (Lattice constant (a.u.))\n")
        f.write(f"  {vec_a[0]:<8.4f}{vec_a[1]:<8.4f}{vec_a[2]:<8.4f}\n")
        f.write(f"  {vec_b[0]:<8.4f}{vec_b[1]:<8.4f}{vec_b[2]:<8.4f}\n")
        f.write(f"  {vec_c[0]:<8.4f}{vec_c[1]:<8.4f}{vec_c[2]:<8.4f}\n")
        f.write(f"   {len(atom_types)}\n")
        for z_number, basis_atoms in atom_types.items():
            symbol = atomic_symbols.get(z_number, 'X')
            f.write(f"{symbol}\n")
            f.write(f"   {len(basis_atoms):<2d}{float(z_number):<8.4f} 0.0000 1.0000\n")
            for atom_coord in basis_atoms:
                f.write(f"  {atom_coord[0]:<8.4f}{atom_coord[1]:<8.4f}{atom_coord[2]:<8.4f}\n")
        f.write("   3\n")
        f.write("  0.6667\n")
        f.write("  10\n") 
    print("Arquivo de input do MSCD gerado com sucesso.")


# --- Bloco Principal MODIFICADO ---
if __name__ == "__main__":
    input_file = 'input_cluster.txt'
    
    # 1. Lê o parâmetro de rede do cabeçalho primeiro
    lattice_multiplier = get_lattice_constant_from_header(input_file)
    
    # 2. Passa o multiplicador para a função que lê as camadas
    parametros = parse_multi_layer_file(input_file, lattice_constant=lattice_multiplier)
    
    # 3. O resto continua como antes
    element_order = get_element_order_from_header(input_file)

    if parametros and element_order:
        generate_mscd_phaseshift_input(parametros, output_filename='cluster.i')
        assemble_mufftin_d(element_order, output_filename='atomic.i')

