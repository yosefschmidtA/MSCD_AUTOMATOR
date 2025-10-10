import subprocess
import re
import os
import glob
import shutil

def split_potential_file(input_filename='mufftin.d', output_prefix='mufftin'):
    """
    Lê o arquivo de potencial gerado 'mufftin.d' e o divide em
    múltiplos arquivos, um para cada bloco de dados de elemento.
    """
    print(f"\n--- Iniciando Etapa 3: Divisão do arquivo {input_filename} ---")
    try:
        with open(input_filename, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"ERRO: Arquivo de potencial '{input_filename}' não encontrado. Etapa de divisão pulada.")
        return

    # CORREÇÃO: Usa uma expressão regular para o delimitador,
    # aceitando qualquer número de elementos (ex: NRR= 3, NRR= 4, etc.).
    # O parêntese cria um grupo de captura, que faz com que re.split mantenha o delimitador.
    delimiter_regex = r'(\s*&NL2 NRR=\s*\d+\s*&END)'
    
    # Divide o conteúdo pelo delimitador, mantendo-o. Remove strings vazias.
    parts = [part for part in re.split(delimiter_regex, content) if part.strip()]

    # A lista 'parts' agora deve ter o formato [DELIM1, DADOS1, DELIM2, DADOS2, ...]
    # O número de elementos é a metade do tamanho da lista.
    num_elements = len(parts) // 2

    if num_elements == 0:
        print(f"Nenhum bloco de dados completo (delimitador + dados) encontrado em '{input_filename}'.")
        # Lógica de fallback para caso o arquivo tenha apenas um bloco sem delimitador inicial (caso raro)
        if len(parts) == 1:
            print("Encontrado apenas 1 bloco de dados. Criando mufftin1.d...")
            output_filename_single = f"{output_prefix}1.d"
            # Neste caso, o delimitador não foi encontrado, então não o adicionamos.
            with open(output_filename_single, 'w') as f_out:
                f_out.write(parts[0])
            print(f"Arquivo '{output_filename_single}' criado com sucesso.")
        return

    print(f"Encontrados {num_elements} blocos de dados. Iniciando a divisão...")
    for i in range(num_elements):
        output_filename = f"{output_prefix}{i+1}.d"
        
        # Pega o delimitador e os dados para o elemento atual
        delimiter = parts[i*2]
        data_block = parts[i*2 + 1]

        # Recria o bloco com seu delimitador original
        block_content = f"{delimiter.strip()}\n{data_block.lstrip()}"
        with open(output_filename, 'w') as f_out:
            f_out.write(block_content)
        with open(output_filename, 'r') as f:
            lines = f.readlines()
        if len(lines) > 1:  # garante que o arquivo tem mais de 1 linha
            with open(output_filename, 'w') as f:
                f.writelines(lines[1:])
            
        print(f"Arquivo '{output_filename}' criado com sucesso.")
    
    print("--- Divisão de arquivos concluída. ---")

def run_poconv_sequence(executable_path='./poconv', input_prefix='mufftin'):
    """
    Encontra todos os arquivos 'mufftin*.d' e executa o programa poconv
    para cada um deles, gerando os arquivos 'ps*.txt'.
    """
    print(f"\n--- Iniciando Etapa 4: Execução do {executable_path} ---")
    if not os.access(executable_path, os.X_OK):
        print(f"ERRO: O arquivo '{executable_path}' não é executável ou não foi encontrado.")
        return

    # CORREÇÃO: Padrão de glob mais específico para pegar apenas os arquivos divididos.
    muffin_files = sorted(glob.glob(f"{input_prefix}[0-9]*.d"))
    if not muffin_files:
        print(f"Nenhum arquivo '{input_prefix}[0-9]*.d' encontrado para processar. Etapa do poconv pulada.")
        return

    print(f"Arquivos encontrados para processar: {', '.join(muffin_files)}")
    for input_file in muffin_files:
        match = re.search(r'(\d+)\.d$', input_file)
        if not match:
            continue
        
        file_number = match.group(1)
        output_file = f"ps{file_number}.txt"
        print(f"\nProcessando '{input_file}' -> '{output_file}'...")

        try:
            # Executa o processo e mostra a saída em tempo real
            process_poconv = subprocess.Popen(
                [executable_path], stdin=subprocess.PIPE, text=True
            )
            poconv_input = f"{input_file}\n{output_file}\n1\n1\n"
            process_poconv.communicate(input=poconv_input, timeout=15)

            if process_poconv.returncode != 0:
                print(f"Ocorreu um erro ao processar '{input_file}'. Verifique a saída do programa acima.")
            else:
                print(f"'{output_file}' criado com sucesso.")
        except subprocess.TimeoutExpired:
            print(f"ERRO: O processo '{executable_path}' demorou muito para responder para o arquivo '{input_file}'.")
    
    print("--- Execução do poconv concluída. ---")

def parse_element_info(filename='input_cluster.txt'):
    """
    Lê as linhas de cabeçalho do 'input_cluster.txt' para extrair
    informações sobre cada elemento (orbital, energia).
    """
    print(f"\n--- Lendo informações de elementos de {filename} ---")
    element_data = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split()
                try:
                    int(parts[0])
                    float(parts[3])
                    is_valid_header = True
                except (ValueError, IndexError):
                    is_valid_header = False

                if len(parts) >= 4 and is_valid_header:
                    info = {
                        'index': int(parts[0]), 'atomic_number': int(parts[1]),
                        'orbital': parts[2], 'energy': parts[3]
                    }
                    element_data.append(info)
                    print(f"  - Encontrado elemento {info['index']}: Z={info['atomic_number']}, Orbital={info['orbital']}, Energia={info['energy']}")
                else:
                    break
        return element_data
    except FileNotFoundError:
        print(f"ERRO: Arquivo de configuração '{filename}' não encontrado.")
        return None

def generate_psrmin_files(element_info_list):
    """
    Gera um arquivo psrmin.txt para cada elemento encontrado.
    """
    print("\n--- Iniciando Etapa 5: Geração dos arquivos psrmin ---")
    if not element_info_list:
        print("Nenhuma informação de elemento para processar. Etapa pulada.")
        return

    psrmin_template = """     821    11    20      datakind beginning-row linenumbers

  ----------------------------------------------------------------
                         David A. Shirley's group
                   Pennsylvania State University (PSU)
               Lawrence Berkeley National Laboratory (LBNL)
     Copyright (c) 1995-1996 DAS group. All rights reserved
  ----------------------------------------------------------------

         input file for phase shift or radial matrix calculation 
'po'    '{ps_input_file}'         input potential data filein
'ps'    '{ps_output_file}'      output phase shift data file
'rm'    '{rm_output_file}'      output radial matrix data file
'ss'    '{orbital}'               subshell and initial state
'sb'    '1'             symbol of atom
'at'    '1'             name of atom

10         1    lnum, output file (0=phase shfit, 1= radial matrix)
5.0         20.00   0.05    kin, kmax, kstep
{energy}                     subshell binding energy
"""
    for info in element_info_list:
        index = info['index']
        output_filename = f"psrmin{index}.txt"
        content = psrmin_template.format(
            ps_input_file=f"ps{index}.txt", ps_output_file=f"ps{index}.1.txt",
            rm_output_file=f"rm{index}.txt", orbital=info['orbital'], energy=info['energy']
        )
        with open(output_filename, 'w') as f:
            f.write(content)
        print(f"Arquivo '{output_filename}' criado com sucesso.")
    print("--- Geração dos arquivos psrmin concluída. ---")

def run_psrm_sequence(executable_path='./psrm.x'):
    """
    Executa o psrm.x para cada elemento, primeiro para phase shift (0)
    e depois para radial matrix (1).
    """
    print(f"\n--- Iniciando Etapa 6: Execução do {executable_path} ---")
    if not os.access(executable_path, os.X_OK):
        print(f"ERRO: O arquivo '{executable_path}' não é executável ou não foi encontrado.")
        return

    # CORREÇÃO: Padrão de glob mais específico para evitar pegar 'psrmin.txt'
    psrmin_files = sorted(glob.glob("psrmin[0-9]*.txt"))
    if not psrmin_files:
        print("Nenhum arquivo 'psrmin[0-9]*.txt' encontrado para processar. Etapa pulada.")
        return

    for psrmin_file in psrmin_files:
        print(f"\nProcessando {psrmin_file}...")
        with open(psrmin_file, 'r') as f:
            original_content = f.read()

        # --- Cálculo do Phase Shift (0) ---
        print("  - Configurando para Phase Shift (0)...")
        content_ps = re.sub(r"^(10\s+)1(\s+lnum.*)$", r"\g<1>0\g<2>", original_content, flags=re.MULTILINE)
        with open('psrmin.txt', 'w') as f:
            f.write(content_ps)
        
        print(f"  - Executando {executable_path} para Phase Shift...")
        subprocess.run([executable_path], text=True)
        print("  - Phase Shift calculado.")

        # --- Cálculo da Radial Matrix (1) ---
        print("  - Configurando para Radial Matrix (1)...")
        content_rm = re.sub(r"^(10\s+)0(\s+lnum.*)$", r"\g<1>1\g<2>", content_ps, flags=re.MULTILINE)
        with open('psrmin.txt', 'w') as f:
            f.write(content_rm)
            
        print(f"  - Executando {executable_path} para Radial Matrix...")
        subprocess.run([executable_path], text=True)
        print("  - Radial Matrix calculada.")

    # Limpa o arquivo temporário no final
    try:
        os.remove('psrmin.txt')
        print("\nArquivo temporário 'psrmin.txt' removido.")
    except FileNotFoundError:
        pass # Não faz nada se o arquivo não existir

    print("--- Execução do psrm.x concluída. ---")


def run_full_sequence(phsh1_exec='./phsh1', poconv_exec='./poconv', psrm_exec='./psrm.x', config_file='input_cluster.txt', potential_file='atomic.i'):
    """
    Orquestra a execução de toda a sequência de automação.
    """
    # ETAPA 1: Cálculo "Bulk" para obter o valor MTZ
    print("--- Iniciando Etapa 1: Cálculo Bulk (opção 0) ---")
    mtz_value = None
    try:
        if not os.access(phsh1_exec, os.X_OK):
            print(f"ERRO: O arquivo '{phsh1_exec}' não é executável.")
            return

        process_bulk = subprocess.Popen(
            [phsh1_exec], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True
        )
        stdout_bulk, stderr_bulk = process_bulk.communicate(input='0\n')
        if process_bulk.returncode != 0:
            print(f"Ocorreu um erro na Etapa 1:\n{stderr_bulk}")
            return
        print(f"Saída do cálculo Bulk:\n{stdout_bulk}")
        lines = stdout_bulk.strip().split('\n')
        if lines:
            last_line = lines[-1].strip()
            try:
                mtz_value = float(last_line)
                print(f"+++ Valor MTZ extraído com sucesso: {mtz_value} +++\n")
            except ValueError:
                pass
        if mtz_value is None:
            print("ERRO: Não foi possível extrair o valor MTZ da saída da Etapa 1.")
            return
    except FileNotFoundError:
        print(f"ERRO: Executável '{phsh1_exec}' não encontrado.")
        return

    # ETAPA 2: Cálculo "Slab" com o valor MTZ
    print("--- Iniciando Etapa 2: Cálculo Slab (opção 1) ---")
    try:
        # A saída desta etapa será mostrada em tempo real
        process_slab = subprocess.Popen(
            [phsh1_exec], stdin=subprocess.PIPE, text=True
        )
        input_slab = f"1\n{mtz_value}\n"
        print(f"Enviando para o programa:\n1\n{mtz_value}")
        process_slab.communicate(input=input_slab)
        if process_slab.returncode != 0:
            print(f"Ocorreu um erro na Etapa 2. Verifique a saída acima.")
            return
        print("\n--- Sequência de cálculo concluída com sucesso! ---")
    except FileNotFoundError:
        print(f"ERRO: Executável '{phsh1_exec}' não encontrado.")
        return

    # ETAPAS 3, 4, 5 e 6
    split_potential_file(input_filename='mufftin.d')
    run_poconv_sequence(executable_path=poconv_exec)
    element_info = parse_element_info(filename=config_file)
    if element_info:
        generate_psrmin_files(element_info)
        run_psrm_sequence(executable_path=psrm_exec)


# --- Bloco Principal ---
if __name__ == "__main__":
    run_full_sequence()

