import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal
from scipy.stats import sem
from spectra_temperature import get_temperature_from_spectra, baseline_correction, find_peaks_position
import os

def calc(NIST_db_name, spectras_name):
    # Read NIST database fo Bi
    NIST_db_folder = 'NIST'
    nist = []
    with open(f'{NIST_db_folder}/{NIST_db_name}.txt', 'r') as f:
        for line in f:
            nist.append(line.split())

    NIST_indexes = [num for num, i in enumerate(nist) if float(i[0]) > 290 and float(i[0]) < 500]
    nist = nist[min(NIST_indexes):max(NIST_indexes)+1]
    NIST_wavelength = [float(x[0]) for x in nist]

    NIST_A = [float(x[4]) for x in nist]
    NIST_E = [float(x[8])*0.000123984 for x in nist]
    NIST_G = [int(int(i)/int(j)*2) + 1 for i,j in [x[14].split('/') for x in nist]]
    NIST_data = [NIST_wavelength, NIST_A, NIST_E, NIST_G]
    #  ----------------------------------------------------------

    # Чтение спектра 
    data = []
    spectras_folder = 'Experimental Data'
    with open(f'{spectras_folder}/{spectras_name}.txt', 'r') as f:
        for line in f:
            line = line.replace(',', '.')
            data.append(line.split())

    # Берем имена колонок фала
    col_names = data.pop(0)

    # Транспонирование прочитанных строк в колонки (теперь каждый список - спектр по длинам волн)
    spectras = np.asarray(data[400:1400], dtype='float')
    spectras = np.transpose(spectras)

    # ось Х в пикселях
    pixel_numbers = spectras[0]
    # ось Х в длинах волн
    wavelength = spectras[1]
    spectras = np.delete(spectras, [0,1,3], 0)

    # ---------------------------
    # Вычитание базовой линии из всех спектров
    baselines = []
    spectras_corrected = []
    for spectra in spectras:
        baseline = baseline_correction(spectra, 50, 0.00001, 200)
        new_spectra = [a-b for a,b in zip(spectra, baseline)]
        baselines.append(baseline)
        spectras_corrected.append(np.asarray(new_spectra))

    # В переменной лежит массив вида (длина волны1, длина волны2, температура рассчитанная по этим линиям)
    all_temperature_data = []
    # РАссчет температур по суперпозиции всех линий из БД NIST
    for i, spectr in enumerate(spectras_corrected):
        result = get_temperature_from_spectra(i, NIST_data, [wavelength, spectr])
        if result != []:
            all_temperature_data.append(result)
    print(all_temperature_data)

    # распределение температур по парам линий с погрешностями
    peak_to_peak_data = {}

    for spec in all_temperature_data:
        for row in spec:
            new_row = sorted(row[0:2])
            new_row.append(row[2])
            if f'{new_row[0]} {new_row[1]}' in peak_to_peak_data:
                peak_to_peak_data[f'{new_row[0]:.1f}-{new_row[1]:.1f}'].append(new_row[2])
            else:
                peak_to_peak_data[f'{new_row[0]:.1f}-{new_row[1]:.1f}'] = [new_row[2]]

    # данные о расчетах температуры по всем комбинациям пиков
    peak_to_peak_temperatures = []
    for key, value in peak_to_peak_data.items():
        peak_to_peak_temperatures.append([key, np.average(value), sem(value)])


    # Список средних температур и ошибкок по каждому спектру
    average_temperatures = []
    average_temperatures_errors = []
    for spectra_temperatures in all_temperature_data:   
        line_by_line_temperatures = []
        for line in spectra_temperatures:
            line_by_line_temperatures.append(line[2])
        average_temperatures.append(np.average(line_by_line_temperatures))
        average_temperatures_errors.append(sem(line_by_line_temperatures))
    final_temperature = np.average(average_temperatures)
    final_error = np.average(average_temperatures_errors)


    # Проверка директории результатов
    if not os.path.exists('results'):
        os.makedirs('results')

    # создание директории отчета
    result_folder = f'results/{NIST_db_name}_{spectras_name}'
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)


    with open(f'{result_folder}/peak_to_peak_temperatures_data.txt', 'w') as f:
        for key, value in peak_to_peak_data.items():
            f.write(f'{key}: {value}\n')

    with open(f'{result_folder}/peak_to_peak_temperatures.txt', 'w') as f:
        for line in peak_to_peak_temperatures:
            f.write(f'{line}\n')

    with open(f'{result_folder}/temperature_data_all.txt', 'w') as f:
        for num, line in enumerate(all_temperature_data):
            f.write(f"----------- spectra N{num}\n")
            for subline in line:
                f.write(f"{subline}\n")

    # Файл с температурами и ошибками для каждого спектра
    with open(f'{result_folder}/temperature_result.txt', 'w') as f:
        for num, temperature_data in enumerate(zip(average_temperatures,average_temperatures_errors)):
            f.write(f"{num} {temperature_data[0]}+-{temperature_data[1]}\n")
        f.write(f"---------------------------------------------------")
        f.write(f"result =  {final_temperature} +- {final_error}\n")


    # -------------------  Графики для каждого спектра -------------------

    # Линия нуля
    zero_y = [0,0]
    zero_x = [280, 540]
    x_tick_label = [300,350,400,450,500,550]
    y_tick_label = [0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000]

    for i in range(len(spectras)):
        plt.title(f'Спектр номер {i}')
        plt.ylabel("Интенсивность")
        plt.xlabel("Длина волны")
        # Линия нуля
        plt.plot(zero_x, zero_y)
        # начальный спектр
        plt.plot(wavelength, spectras[i])
        # базовая линия
        plt.plot(wavelength, baselines[i])
        # спектр с вычтенной базовой линией
        plt.plot(wavelength, spectras_corrected[i])
        # Отмечаем пики из NIST исользуемые в рассчетах
        used_peaks_positions = find_peaks_position(NIST_wavelength, [wavelength, spectras_corrected[i]])
        plt.plot(NIST_wavelength , spectras_corrected[i][used_peaks_positions], "X")
        # Убираем подписи осей (мешают пока)
        plt.xticks(x_tick_label) 
        plt.yticks(y_tick_label) 
        plt.savefig(f'{result_folder}/spectra_plot{i}.png')
        plt.clf()

        with open(f'{result_folder}/spectra_data{i}.txt', 'w') as f:
            for text_spectra in zip(wavelength, spectras[i], baselines[i], spectras_corrected[i]):
                f.write(f"{text_spectra[0]} {text_spectra[1]} {text_spectra[2]} {text_spectra[3]}\n")

    # Графики зависимости температуры от номера спектра, средняя температура
    index = range(0,len(average_temperatures))
    plt.bar(index, average_temperatures, yerr = average_temperatures_errors, error_kw={'ecolor':'0.1','capsize':6},alpha=0.7, label=f'avg T = {final_temperature} +- {final_error}')
    plt.legend(loc=2)
    plt.title('Температура плазмы от номера спектра')
    plt.ylabel("Температура плазмы")
    plt.xlabel("Номер Спектра")
    plt.savefig(f'{result_folder}/average_temperature.png')
    plt.clf()

    # Графики зависимости температуры от выбора линий в спектре
    lines_used = [x[0] for x in peak_to_peak_temperatures]
    lines_temperature = [x[1] for x in peak_to_peak_temperatures]
    lines_error = [x[2] for x in peak_to_peak_temperatures]
    plt.bar(lines_used, lines_temperature, yerr = lines_error, error_kw={'ecolor':'0.1','capsize':6},alpha=0.7, label=f'avg T = {np.average(lines_temperature)} +- {sem(lines_temperature)}')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.legend(loc=2)
    plt.title('Температура плазмы от выбора линии в спектре')
    plt.ylabel("Температура плазмы")
    plt.xlabel("Пара линий")
    plt.savefig(f'{result_folder}/average_temperature_by_lines.png')
    plt.clf()

    return spectras_name, final_temperature, final_error


if __name__ == '__main__':
        
    np.seterr(all="ignore")
    energies = ['8', '8-5', '9', '9-5', '10', '10-5', '11', '11-5', '12']
    #energies = ['12']
    res = {}
    for i in energies:
        spectras_name, final_temperature, final_error = calc('Bi', i)
        res[spectras_name] = [final_temperature, final_error]

    x =[]
    y = []
    err = []
    for name, value in res.items():
        print(f'{name}: {value}')
        x.append(name)
        y.append(value[0])
        err.append(value[1]/2)


    plt.ylim([6000, 8000])
    plt.title(f'Спектр номер {i}')
    plt.ylabel("Интенсивность")
    plt.xlabel("Длина волны")
    # Линия нуля
    print(x, y, err)
    plt.errorbar(x, y, yerr=err, fmt='o')
    # начальный спектр

    plt.savefig(f'results/res_long.png')
    plt.clf()