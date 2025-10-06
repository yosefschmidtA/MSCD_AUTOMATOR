import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
def calcular_r_factor(df):
    resultados = []
    angulos = []
    resultado_final4 = 0
    n_total = 0

    angulos_unicos = df['Theta'].unique()

    # Acumuladores para o cálculo do Rtotal
    todas_intensidades_exp = []
    todas_intensidades_calc = []

    for theta in angulos_unicos:
        subset = df[df['Theta'] == theta]
        phi = subset['Phi'].values
        intensidade_calculada = np.round(subset['intensitycal'].values, 10)
        intensidade_experimental = np.round(subset['intensityexp'].values, 10)

        # Acumular todas as intensidades para o cálculo do Rtotal
        todas_intensidades_exp.extend(intensidade_experimental)
        todas_intensidades_calc.extend(intensidade_calculada)

        # Normalização por theta
        somatorio1 = np.round(np.sum(np.abs(intensidade_experimental)),10)
        somatorio2 = np.round(np.sum(np.abs(intensidade_calculada)),10)

        if somatorio1 == 0 or somatorio2 == 0:
            continue

        nova_coluna_exp = np.round((intensidade_experimental / somatorio1),10)
        nova_coluna_calc = np.round((intensidade_calculada / somatorio2),10)

        resultado_final1 = np.round(np.sum((nova_coluna_exp - nova_coluna_calc) ** 2),10)
        resultado_final2 = np.round(np.sum(nova_coluna_exp ** 2 + nova_coluna_calc ** 2),10)

        resultado_final3 = np.round((resultado_final1 / resultado_final2),10)
        resultado_final4 += resultado_final3

        resultados.append(resultado_final3)
        angulos.append(theta)

        n_total += 1

    # Cálculo do R-factor médio
    r_factor_medio = np.round((resultado_final4 / n_total),10) if n_total > 0 else None

    # Cálculo do Rtotal usando todas as intensidades acumuladas
    todas_intensidades_exp = np.round(np.array(todas_intensidades_exp),10)
    todas_intensidades_calc = np.round(np.array(todas_intensidades_calc),10)

    somatorio1_total = np.round(np.sum(np.abs(todas_intensidades_exp)),10)
    somatorio2_total = np.round(np.sum(np.abs(todas_intensidades_calc)),10)

    if somatorio1_total > 0 and somatorio2_total > 0:
        nova_coluna_exp_total = np.round((todas_intensidades_exp / somatorio1_total),10)
        nova_coluna_calc_total = np.round((todas_intensidades_calc / somatorio2_total),10)

        resultado_final1_total = np.round(np.sum((nova_coluna_exp_total - nova_coluna_calc_total) ** 2),10)
        resultado_final2_total = np.round(np.sum(nova_coluna_exp_total ** 2 + nova_coluna_calc_total ** 2),10)

        r_factor_total = np.round((resultado_final1_total / resultado_final2_total),10)
    else:
        r_factor_total = None

    return resultados, angulos, r_factor_medio, r_factor_total
def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    theta_value = None
    lista_temporaria = []

    for i in range(26, len(lines)):
        line = lines[i].strip()

        if "fitted parameters (" in line:
            break

        if line:
            parts = line.split()

            if len(parts) == 7:
                # Encontrou um novo Theta, processa os dados acumulados
                if lista_temporaria:
                    intensi2_values = [item[1] for item in lista_temporaria]
                    media_intensi2 = sum(intensi2_values) / len(intensi2_values)

                    for phi, intensi1, intensi2, intensityexp in lista_temporaria:
                        intensitycal = (intensi1 - media_intensi2) / media_intensi2
                        data.append([phi, intensitycal, theta_value, intensityexp, True])

                    lista_temporaria = []

                theta_value = float(parts[3])

            elif len(parts) == 5 and theta_value is not None:
                phi = float(parts[0])
                intensi1 = float(parts[1])
                intensi2 = float(parts[2])
                intensityexp = float(parts[4])
                lista_temporaria.append((phi, intensi1, intensi2, intensityexp))

    # Depois do loop, processa o último conjunto também
    if lista_temporaria:
        intensi2_values = [item[1] for item in lista_temporaria]
        media_intensi2 = sum(intensi2_values) / len(intensi2_values)

        for phi, intensi1, intensi2, intensityexp in lista_temporaria:
            intensitycal = (intensi1 - media_intensi2) / media_intensi2
            data.append([phi, intensitycal, theta_value, intensityexp, True])

    df = pd.DataFrame(data, columns=['Phi', 'intensitycal', 'Theta', 'intensityexp', 'IsOriginal'])
    resultados, angulos, r_factor_medio, r_factor_total = calcular_r_factor(df)
    for theta, r_factor in zip(angulos, resultados):
        print(f'Theta: {theta}, R-factor: {r_factor}')

    print(f'R-factor médio: {r_factor_medio}')
    print(f'R-factor total: {r_factor_total}')
    # Verificar o intervalo de Phi




    if not np.isclose(df['Phi'].min(), 0):
        df_0 = df[df['Phi'] == df['Phi'].min()].copy()
        df_0['Phi'] = 0
        df = pd.concat([df, df_0], ignore_index=True)

    phi_min = df['Phi'].min()
    phi_max = df['Phi'].max()
    phi_interval = phi_max - phi_min

    if phi_interval < 360 and df['Phi'].max() < 360:
        df_360 = df[df['Phi'] == 0].copy()
        df_360['Phi'] = 360
        df = pd.concat([df, df_360], ignore_index=True)

    if phi_interval == 120:
        # Replicação dos dados para cobrir 360 graus
        first_values = df.groupby('Theta').first().reset_index()
        df = df.groupby('Theta', group_keys=False, as_index=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original
        # Adicionar os novos valores ao DataFrame
        df = pd.concat([df, last_values], ignore_index=True)

        df_0_120 = df.copy()
        df_0_120['Phi'] = 120 + df_0_120['Phi']
        df_0_120['isOriginal'] = False

        df_240_360 = df.copy()
        df_240_360['Phi'] = 240 + df_240_360['Phi']

        df = pd.concat([df, df_0_120, df_240_360]).reset_index(drop=True)

    if phi_interval == 90:
        # Replicação dos dados para cobrir 360 graus
        first_values = df.groupby('Theta').first().reset_index()
        df = df.groupby('Theta', group_keys=False, as_index=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original

        # Adicionar os novos valores ao DataFrame
        df = pd.concat([df, last_values], ignore_index=True)
        df_0_90 = df.copy()
        df_0_90['Phi'] = 90 + df_0_90['Phi']

        df_90_180 = df.copy()
        df_90_180['Phi'] = 180 + df_90_180['Phi']

        df_180_270 = pd.concat([df,df_0_90]).reset_index(drop=True)
        df_180_270['Phi'] = 180 + df_180_270['Phi']

        df = pd.concat([df, df_0_90, df_180_270]).reset_index(drop=True)

    return df, r_factor_total  # Retorna r_factor_total

def interpolate_data(df, resolution=1000):
    phi = np.radians(df['Phi'])
    theta = np.radians(df['Theta'])
    intensity = df['intensitycal']

    phi_grid = np.linspace(np.min(phi), np.max(phi), resolution)
    theta_grid = np.linspace(np.min(theta), np.max(theta), resolution)

    phi_grid, theta_grid = np.meshgrid(phi_grid, theta_grid)

    intensity_grid = griddata((phi, theta), intensity, (phi_grid, theta_grid), method='cubic')

    return phi_grid, theta_grid, intensity_grid


def plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=None, save_path=None):
    import matplotlib.pyplot as plt
    import numpy as np

    # 1. Manter apenas a banda real esquerda (90°–270°)
    df_left = df[(df['Phi'] >= 90) & (df['Phi'] <= 270)].copy()

    # 2. Criar a banda direita por espelhamento em 180°
    df_right = df_left.copy()
    df_right['Phi'] = (df_right['Phi'] + 180) % 360

    # 3. Eliminar os dados reais da direita (para evitar sobreposição)
    df_clean = df[(df['Phi'] >= 90) & (df['Phi'] <= 270)].copy()

    # 4. Montar o DataFrame final com simetria
    df_plot = pd.concat([df_clean, df_right], ignore_index=True)

    def rotate_phi_for_plot(df_plot, rotation_angle):
        """
        Mesma rotação, mas garante que Phi = 0 esteja presente **apenas** para plotagem.
        """
        
        df_plot['Phi'] = (df_plot['Phi'] + rotation_angle) % 360

        if not np.isclose(df_plot['Phi'].min(), 0):
            df_0 = df_plot[df_plot['Phi'] == df_plot['Phi'].min()].copy()
            df_0['Phi'] = 0
            df_plot = pd.concat([df_plot, df_0], ignore_index=True)

        return df_plot
    df_plot = rotate_phi_for_plot(df_plot, 0)

    # Interpolação e gráfico polar
    plt.ion()
    phi_grid, theta_grid, intensity_grid = interpolate_data(df_plot, resolution)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 8), dpi=100)
    c = ax.pcolormesh(phi_grid, theta_grid, intensity_grid, shading='gouraud', cmap='afmhot')

    max_theta = df['Theta'].max()
    ax.set_ylim(0, np.radians(max_theta))
    theta_ticks = np.linspace(0, max_theta, num=6)
    ax.set_yticks(np.radians(theta_ticks))
    ax.set_yticklabels([f'{int(t)}°' for t in theta_ticks])

    cbar = fig.colorbar(c, ax=ax, label='', pad=0.08)
    cbar.set_label('', fontsize=22, fontweight='bold')
    cbar.ax.yaxis.set_tick_params(labelsize=30)
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')

    if my_variable is not None:
        fig.text(0.87, 0.03, f'R-factor: {my_variable}', fontsize=26, color='black', ha='right', va='bottom', fontweight='bold')

    fig.text(0.94, 0.9, "Anisotropy", fontsize=34, color='black', ha='right', va='bottom', fontweight='bold')

    phi_ticks = np.linspace(0, 2 * np.pi, num=9)[:-1]
    phi_labels = [f'{int(np.degrees(tick))}°' for tick in phi_ticks]
    ax.set_xticks(phi_ticks)
    ax.set_xticklabels(phi_labels, fontsize=26, fontweight='bold')

    pad_values = [1, -1, 3, 0, -7, -6, -1, -6]
    for label, pad in zip(ax.get_xticklabels(), pad_values):
        label.set_y(label.get_position()[1] + pad * 0.01)

    ax.tick_params(pad=8)
    plt.yticks(fontsize=0, fontweight='bold')
    plt.draw()

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')

    plt.pause(600)


# Caminho do arquivo
file_path = './arquivos/saida.out'
save_path = 'grafico_p.png'

df, r_factor_total = process_file(file_path)
#r_factor_total=0.276
my_variable = "{:.3f}".format(r_factor_total)
plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=my_variable, save_path=save_path)
