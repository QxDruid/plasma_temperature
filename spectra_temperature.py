from cmath import nan
import numpy as np
import math
from scipy import signal, sparse, linalg


# Поиск базовой линии методом Asymmetric Least Squares Smoothing
def baseline_correction(y, lam, p, iter=10):
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(iter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sparse.linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

# Поиск реальных положений линий в спектре по базе NIST. параметр spectra = кортеж из длины волны и интенсивности ( то есть х и у оси спектра)
def find_peaks_position(NIST_wavelength, spectra):
    wavelength, spectra_intensity = spectra
    peak_positions = []
    for peak in NIST_wavelength:
        x = [abs(i - peak) for i in wavelength]
        peak_positions.append(np.argmin(x))

    # Уточнение максимумов линий (реальные максимумы могут отличаться на пару пикселей от данных NIST)
    peak_positions_real = []
    const_sdvig = 3
    for peak in peak_positions:
        peak_positions_real.append(np.argmax(spectra_intensity[peak-const_sdvig:peak+const_sdvig]) + peak - const_sdvig)

    return peak_positions_real

# Функция рассчета температуры плазмы по суперпозиции всех пиков вещества из базы данных NIST
def get_temperature_from_spectra(i, NIST_data, spectra):
    wavelength, spectra_intensity = spectra
    spectra_intensity = np.asarray(spectra_intensity)
    wavelength = np.asarray(wavelength)
    NIST_wavelength, NIST_A, NIST_E, NIST_G = NIST_data
    
    # Поиск линий висмута в реальном спектре с S100
    peak_positions_real = find_peaks_position(NIST_wavelength, spectra)
    # Интенсивность найденных линий
    lines_intensity = spectra_intensity[peak_positions_real]
    # Собираем все параметры для каждой линии в один лист
    line_nums = [num for num in range(len(NIST_wavelength))]
    parameter_list = list(zip(line_nums, NIST_wavelength, lines_intensity, NIST_A, NIST_E, NIST_G))

    # убираем слишком маленькие линии
    avg_intensity = np.average(lines_intensity)
    parameter_list = [x for x in parameter_list if x[2] > avg_intensity/5]

    k = 0.00008617333262
    temperature = []
    for num, line1 in enumerate(parameter_list):
        parameter_list.pop(num)
        for line2 in parameter_list:
            if line1[0] != line2[0] and line1[2] > 0 and line2[2] > 0:
                T = ((line1[4] - line2[4])/k) * 1/(math.log((line1[2]*line2[3]*line2[5]*line1[1]) / (line2[2]*line1[3]*line1[5]*line2[1])))
                print(f't= {T}') 
                if abs(T) > 0:
                    temperature.append((line1[1], line2[1], abs(T)))
                else:
                    pass
                    #print(f' false temp {(line1[1], line2[1], abs(T))}')
    print(temperature)
    return temperature
