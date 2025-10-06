import numpy as np
import math
import re

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
    """Lê um arquivo com múltiplas layers separadas por '///'."""
    print(f"Lendo o arquivo de configuração do cluster: {filename}")
    try:
        with open(filename, 'r') as f:
            content = f.read()

        # Divide o arquivo em blocos de layers e descarta o primeiro bloco (cabeçalho)
        layer_blocks = re.split(r'/{3,}', content)
        data_blocks = layer_blocks[1:]  # Pula o primeiro bloco de cabeçalho

        all_layers_params = []
        for block in data_blocks:
            if not block.strip():
                continue

            params = {}
            lines = [line for line in block.strip().split('\n') if line.strip()]

            # Pega a primeira linha do bloco de dados da camada
            parts = lines[0].split()
            params['atomic_number'] = int(parts[0])
            # A verificação de '1' ou '0' para 'is_emitter_layer' deve ser feita no segundo elemento, não no primeiro
            # A linha original `26 0 (1 ou 0, emitindo ou não emitindo)` pode confundir o parser, que leria `0`
            # O ideal é que a linha tenha apenas `26 0` ou `26 1`
            # Assumindo que o segundo elemento é o flag do emissor
            params['is_emitter_layer'] = (int(parts[1]) == 1)

            for i, line in enumerate(lines):
                if 'unita' in line:
                    parts = lines[i + 1].split()
                    params['unita_len'] = float(parts[0])
                    params['unita_ang_deg'] = float(parts[1])
                elif 'unitb' in line:
                    parts = lines[i + 1].split()
                    params['unitb_len'] = float(parts[0])
                    params['unitb_ang_deg'] = float(parts[1])
                elif 'origin' in line:
                    parts = lines[i + 1].split()
                    params['origin_len'] = float(parts[0])
                    params['origin_ang_deg'] = float(parts[1])
                elif 'interlayer' in line:
                    parts = lines[i + 1].split()
                    params['interlayer_spacing'] = float(parts[0])

            all_layers_params.append(params)

        print(f"Encontradas {len(all_layers_params)} layers para processar.")
        return all_layers_params

    except FileNotFoundError:
        print(f"ERRO: Arquivo de input '{filename}' não foi encontrado!")
        return None
    except Exception as e:
        print(f"ERRO ao ler o arquivo de input: {e}")
        return None

# Mantenha o restante do código igual
def create_3d_cluster_xyz(all_params, radius_angstroms=15.0, output_filename='cluster_3d.xyz'):
    # ... (código da sua função create_3d_cluster_xyz aqui)
    all_atom_coords = []
    accumulated_z = 0.0

    for layer_idx, params in enumerate(all_params):
        print(f"\nProcessando Layer {layer_idx + 1}...")

        atom_symbol = atomic_symbols.get(params['atomic_number'], 'X')
        emitter_symbol = 'He'

        accumulated_z -= params['interlayer_spacing']
        current_z = accumulated_z

        print(f"  - Interlayer Spacing: {params['interlayer_spacing']:.4f} Å")
        print(f"  - Posição Z da camada: {current_z:.4f} Å")

        origin_ang_rad = np.radians(params['origin_ang_deg'])
        unita_ang_rad = np.radians(params['unita_ang_deg'])
        unitb_ang_rad = np.radians(params['unitb_ang_deg'])

        origin_vec = np.array(
            [params['origin_len'] * np.cos(origin_ang_rad), params['origin_len'] * np.sin(origin_ang_rad)])
        unita_vec = np.array([params['unita_len'] * np.cos(unita_ang_rad), params['unita_len'] * np.sin(unita_ang_rad)])
        unitb_vec = np.array([params['unitb_len'] * np.cos(unitb_ang_rad), params['unitb_len'] * np.sin(unitb_ang_rad)])

        layer_2d_coords = []
        max_rep_a = int(math.ceil(radius_angstroms / params['unita_len'])) + 2 if params['unita_len'] > 1e-6 else 0
        max_rep_b = int(math.ceil(radius_angstroms / params['unitb_len'])) + 2 if params['unitb_len'] > 1e-6 else 0

        for m in range(-max_rep_a, max_rep_a + 1):
            for n in range(-max_rep_b, max_rep_b + 1):
                point = origin_vec + m * unita_vec + n * unitb_vec
                if np.linalg.norm(point) <= radius_angstroms:
                    layer_2d_coords.append(point)

        unique_coords_set = {tuple(np.round(coord, 6)) for coord in layer_2d_coords}

        central_atom_coord = None
        if params['is_emitter_layer']:
            min_dist_to_origin = float('inf')
            for coord_tuple in unique_coords_set:
                dist = np.linalg.norm(np.array(coord_tuple) - origin_vec)
                if dist < min_dist_to_origin:
                    min_dist_to_origin = dist
                    central_atom_coord = coord_tuple
            print(
                f"  - Esta é uma camada emissora. Átomo central em ({central_atom_coord[0]:.4f}, {central_atom_coord[1]:.4f}) será marcado.")

        for coord_tuple in unique_coords_set:
            symbol = emitter_symbol if params['is_emitter_layer'] and np.allclose(np.array(coord_tuple), np.array(
                central_atom_coord)) else atom_symbol
            all_atom_coords.append({'symbol': symbol, 'coords': np.array([coord_tuple[0], coord_tuple[1], current_z])})

    print(f"\nGerando arquivo final '{output_filename}' com {len(all_atom_coords)} átomos...")
    with open(output_filename, 'w') as f_out:
        f_out.write(f"{len(all_atom_coords)}\n")
        f_out.write(f"Cluster 3D gerado a partir de {len(all_params)} layers\n")

        sorted_coords = sorted(all_atom_coords, key=lambda d: (d['coords'][2], d['coords'][1], d['coords'][0]))

        for atom_data in sorted_coords:
            c = atom_data['coords']
            f_out.write(f"{atom_data['symbol']:<4s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}\n")

    print("Arquivo .xyz 3D gerado com sucesso!")

# --- Bloco Principal de Execução ---
if __name__ == "__main__":
    parametros_todas_as_camadas = parse_multi_layer_file('input_cluster.txt')

    if parametros_todas_as_camadas:
        create_3d_cluster_xyz(
            parametros_todas_as_camadas,
            radius_angstroms=8.0,
            output_filename='cluster.xyz'
        )
