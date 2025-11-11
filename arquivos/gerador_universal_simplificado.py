import numpy as np
import math
import re

# --- Constantes de Conversão ---
A_TO_BOHR = 1.88973

# Dicionário para converter NÚMERO ATÔMICO para SÍMBOLO
atomic_symbols = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    74: 'W'
}

def parse_multi_layer_file(filename='input_cluster.txt'):
    """Lê um arquivo com múltiplas layers, pulando as 4 primeiras linhas de cabeçalho."""
    print(f"Lendo o arquivo de configuração do cluster: {filename}")
    try:
        with open(filename, 'r') as f:
            all_lines = f.readlines()

        # ALTERAÇÃO: Pula as 4 primeiras linhas que são para o psrm.x e processa o resto
        content = "".join(all_lines[4:])

        layer_blocks = re.split(r'/{3,}', content)
        all_layers_params = []
        for block in layer_blocks:
            if not block.strip(): continue
            params = {}
            lines = [line for line in block.strip().split('\n') if line.strip()]
            parts = lines[0].split()
            params['atomic_number'] = int(parts[0])
            params['is_emitter'] = (int(parts[1]) == 1)
            for i, line in enumerate(lines):
                if 'unita' in line:
                    params.update(
                        {'unita_len': float(lines[i + 1].split()[0]), 'unita_ang_deg': float(lines[i + 1].split()[1])})
                elif 'unitb' in line:
                    params.update(
                        {'unitb_len': float(lines[i + 1].split()[0]), 'unitb_ang_deg': float(lines[i + 1].split()[1])})
                elif 'origin' in line:
                    params.update({'origin_len': float(lines[i + 1].split()[0]),
                                   'origin_ang_deg': float(lines[i + 1].split()[1])})
                elif 'interlayer' in line:
                    params.update({'interlayer_spacing': float(lines[i + 1].split()[0])})
            all_layers_params.append(params)
        print(f"Encontradas {len(all_layers_params)} layers para processar.")
        return all_layers_params
    except Exception as e:
        print(f"ERRO ao ler o arquivo de input: {e}")
        return None


# Cole aqui a sua função create_cluster_xyz completa da resposta anterior
def create_cluster_xyz(all_params, radius_angstroms=15.0, output_filename='cluster_3d.xyz'):
    """Gera um arquivo .xyz 3D explícito para visualização."""
    print(f"\nGerando arquivo XYZ: {output_filename}")
    # ... (cole a função completa aqui) ...
    print("Arquivo .xyz gerado com sucesso.")


def generate_mscd_phaseshift_input(all_params, output_filename='mscd_phaseshift.inp'):
    """Gera o arquivo de input para o MSCD com vetores de rede totalmente adaptativos."""
    print(f"\nGerando arquivo de input para o MSCD: {output_filename}")

    # --- 1. Parâmetros da Rede (LÓGICA TOTALMENTE ADAPTATIVA) ---
    first_layer = all_params[0]

    # Parâmetro de rede base (a_lat) é o menor comprimento de vetor diferente de zero
    lattice_constant_A = min(l['unita_len'] for l in all_params if l['unita_len'] > 0)
    lattice_constant_au = lattice_constant_A * A_TO_BOHR

    # Calcula os vetores cartesianos para a primeira camada
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

    # Normaliza os vetores cartesianos pelo parâmetro de rede para obter os vetores do MSCD
    vec_a = np.array([unita_cartesian[0] / lattice_constant_A, unita_cartesian[1] / lattice_constant_A, 0.0])
    vec_b = np.array([unitb_cartesian[0] / lattice_constant_A, unitb_cartesian[1] / lattice_constant_A, 0.0])
    vec_c = np.array([0.0, 0.0, 3.0])  # Mantém fixo como solicitado

    print("Parâmetros de rede calculados:")
    print(f"  - Parâmetro de rede base (a_lat): {lattice_constant_A:.4f} Å")
    print(f"  - Vetor A (normalizado): ({vec_a[0]:.4f}, {vec_a[1]:.4f}, {vec_a[2]:.4f})")
    print(f"  - Vetor B (normalizado): ({vec_b[0]:.4f}, {vec_b[1]:.4f}, {vec_b[2]:.4f})")

    # --- 2. Agrupar Átomos Base por Tipo ---
    atom_types = {}
    accumulated_z = 0.0
    for params in all_params:
        z_number = params['atomic_number']
        if z_number not in atom_types:
            atom_types[z_number] = []
        origin_ang_rad = np.radians(params['origin_ang_deg'])
        origin_x = params['origin_len'] * np.cos(origin_ang_rad)
        origin_y = params['origin_len'] * np.sin(origin_ang_rad)
        accumulated_z -= params['interlayer_spacing']
        basis_atom = np.array([
            origin_x / lattice_constant_A,
            origin_y / lattice_constant_A,
            accumulated_z / lattice_constant_A
        ])
        atom_types[z_number].append(basis_atom)

    # --- 3. Escrever o Arquivo de Saída ---
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


# --- Bloco Principal ---
if __name__ == "__main__":
    # Cole a função create_cluster_xyz aqui
    input_file = 'input_cluster.txt'
    parametros = parse_multi_layer_file(input_file)
    if parametros:
        generate_mscd_phaseshift_input(parametros, output_filename='cluster.i')

