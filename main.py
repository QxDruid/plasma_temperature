import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal
from scipy.stats import sem
from spectra_temperature import get_temperature_from_spectra, baseline_correction

# Read NIST database fo Bi
filename = 'Bi_NIST.txt'
nist = []
with open(filename, 'r') as f:
    for line in f:
        nist.append(line.split())

NIST_indexes = [num for num, i in enumerate(nist) if float(i[0]) > 290 and float(i[0]) < 500]
nist = nist[min(NIST_indexes):max(NIST_indexes)+1]
NIST_wavelength = [float(x[0]) for x in nist]
print('-'*10)
print(NIST_wavelength)
print('-'*10)
NIST_A = [float(x[4]) for x in nist]
NIST_E = [float(x[8])*0.000123984 for x in nist]
NIST_G = [int(int(i)/int(j)*2) + 1 for i,j in [x[14].split('/') for x in nist]]
NIST_data = [NIST_wavelength, NIST_A, NIST_E, NIST_G]
#  ----------------------------------------------------------

# Чтение спектра 
data = []
filename = '7_5.txt'
with open(filename, 'r') as f:
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
baselines = []
spectras_corrected = []
for spectra in spectras:
    baseline = baseline_correction(spectra, 50, 0.00001, 200)
    new_spectra = [a-b for a,b in zip(spectra, baseline)]
    baselines.append(baseline)
    spectras_corrected.append(new_spectra)

# все полученные данные со всех спектров
all_temperature_data = []
for i, spectr in enumerate(spectras_corrected):
    all_temperature_data.append(get_temperature_from_spectra(i, NIST_data, [wavelength, spectr]))


# распределение температур по парам линий с погрешностями
peak_to_peak_data = {}
with open('re.txt', 'w') as f:
    for spec in all_temperature_data:
        for row in spec:
            new_row = sorted(row[0:2])
            new_row.append(row[2])
            if f'{new_row[0]} {new_row[1]}' in peak_to_peak_data:
                peak_to_peak_data[f'{new_row[0]} {new_row[1]}'].append(new_row[2])
            else:
                peak_to_peak_data[f'{new_row[0]} {new_row[1]}'] = [new_row[2]]

peak_to_peak_temperatures = []
for key, value in peak_to_peak_data.items():
    peak_to_peak_temperatures.append([key, np.average(value), sem(value)])

for line in peak_to_peak_temperatures:
    print(line)



with open('res.txt', 'w') as f:
    for num, line in enumerate(all_temperature_data):
        f.write(f"----------- spectra N{num}\n")
        for subline in line:
            f.write(f"{subline}\n")

# Список средних температур по каждому спектру
average_temperatures = []

for spectra_temperatures in all_temperature_data:   
    temps = []
    for line in spectra_temperatures:
        temps.append(line[2])
    average_temperatures.append(np.average(temps))
energy_temperature = np.average(average_temperatures)
error = sem(average_temperatures)
print('<---------------------------- Final results <----------------------------')
print("Average temps for all spectra")
print(average_temperatures)
print('Average temperature for energy')
print(f'{energy_temperature} +- {error}')
# <----------------------    Графики   ----------------------------------->
'''
# Линия нуля
zero_y = [0,0]
zero_x = [280, 540]
plt.plot(zero_x, zero_y)
# начальный спектр
plt.plot(wavelength, spectras[0])
# базовая линия
plt.plot(wavelength, spline)
# спектр с вычтенной базовой линией
plt.plot(wavelength, spectras_processed)
# пики элемента из NIST
plt.plot(wavelength[peak_positions_real], spectras_processed[peak_positions_real], "x")
# пики испольхуемые в рассчетах
plt.plot(peak_used_positions , peak_used_intensity, "o")
# Убираем подписи осей (мешают пока)
plt.xticks([]) 
plt.yticks([]) 
plt.show()

'''




