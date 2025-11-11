import os

def process_input_cluster(input_file, output_file):
    """
    Processa o arquivo input_cluster.txt para reformatar os dados de acordo
    com o formato de saída desejado e adiciona o resultado a um arquivo.

    Args:
        input_file (str): O nome do arquivo de entrada (input_cluster.txt).
        output_file (str): O nome do arquivo de saída.
    """
    if not os.path.exists(input_file):
        return f"ERRO: O ARQUIVO '{input_file}' NAO FOI ENCONTRADO."

    atomic_num_to_kind = {}
    output_lines = []
    current_layer = 1
    
    with open(input_file, 'r') as infile:
        # 1. Lê as 4 primeiras linhas para o mapeamento
        for _ in range(4):
            line = infile.readline().strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    kind = int(parts[0])
                    atomic_num = int(parts[1])
                    atomic_num_to_kind[atomic_num] = kind
        
        # 2. Pula as próximas 12 linhas (bloco de parâmetros)
        for _ in range(12):
            infile.readline()
        
        # 3. Processa cada bloco de dados separado por '///'
        while True:
            # Lê o cabeçalho do bloco (e descarta o restante da linha)
            header_line = infile.readline()
            if not header_line:
                break
            
            parts = header_line.split()
            if not parts or not all(p.isdigit() for p in parts[:2]):
                continue
            
            atomic_num = int(parts[0])
            emitter = int(parts[1])
            
            # Lê as próximas linhas para extrair os dados
            infile.readline()  # Pula "unita(len ang)"
            unita_values = infile.readline().strip()
            infile.readline()  # Pula "unitb(len ang)"
            unitb_values = infile.readline().strip()
            infile.readline()  # Pula "origin(len,ang) (in Angstrons)"
            origin_values = infile.readline().strip()
            infile.readline()  # Pula "interlayer spacing (A)"
            interlayer_spacing = infile.readline().strip()
            infile.readline()  # Pula "///"

            # 4. Encontra o 'kind' e formata a saída
            kind = atomic_num_to_kind.get(atomic_num, atomic_num)
            
            # Constrói o novo cabeçalho e as linhas de dados no formato de saída
            output_lines.append(f"{current_layer}\t{kind}\t{emitter}\t0\t\t\tlayer, kind, emitter, lineatom\n")
            output_lines.append(f"0\t0\t0\t0\t\t\tlatoms(xa,xb,ya,yb)\n")
            output_lines.append(f" {unita_values}\t\t\tunita(len ang) (fcc (111) structure)\n")
            output_lines.append(f" {unitb_values}\t\t\tunitb(len ang) (in unit of lattice)\n")
            output_lines.append(f" {origin_values}\t\t\torigin(len,ang) (in unit of lattice)\n")
            output_lines.append(f" {interlayer_spacing}\t\t\t\t\t\tinterlayer spacing (unit lattice)\n")
            output_lines.append(f"0.0\t0.0\t0.0\t\t\t\tfit try for spacing, length and units\n")
            output_lines.append("\n") # Adiciona uma quebra de linha entre os blocos
            
            current_layer += 1
            
            # Pula linhas em branco até o próximo bloco
            while True:
                pos = infile.tell()
                next_line = infile.readline()
                if not next_line:
                    break
                if next_line.strip() or next_line.strip() == '///':
                    infile.seek(pos)
                    break
        
        # 5. Adiciona os dados ao arquivo de saída, pulando uma linha antes
        with open(output_file, 'a') as outfile:
            if output_lines:  # Verifica se há conteúdo para escrever
                outfile.write("\n")
                outfile.writelines(output_lines)
            
    return f"DADOS PROCESSADOS E ADICIONADOS A '{output_file}'."

# --- Configuração e Execução ---
INPUT_FILE = "input_cluster.txt"
OUTPUT_FILE = "output_header.txt"

result_message = process_input_cluster(INPUT_FILE, OUTPUT_FILE)
print(result_message)
