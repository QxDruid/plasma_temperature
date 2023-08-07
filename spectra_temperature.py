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


# Функция рассчета температуры плазмы по суперпозиции всех пиков вещества из базы данных NIST
def get_temperature_from_spectra(i, NIST_data, spectra):
    wavelength, spectra = spectra
    spectra = np.asarray(spectra)
    wavelength = np.asarray(wavelength)
    NIST_wavelength, NIST_A, NIST_E, NIST_G = NIST_data
    
    '''
    with open(f'res{i}.txt', 'w') as f:
        for line in spectras_processed:
            f.write(f"{line}\n")
    '''

    # Поиск линий висмута в реальном спектре с S100
    peak_positions = []
    for peak in NIST_wavelength:
        x = [abs(i - peak) for i in wavelength]
        peak_positions.append(np.argmin(x))

    # Уточнение максимумов линий (реальные максимумы могут отличаться на пару пикселей от данных NIST)
    peak_positions_real = []
    const_sdvig = 3
    for peak in peak_positions:
        peak_positions_real.append(np.argmax(spectra[peak-const_sdvig:peak+const_sdvig]) + peak - const_sdvig)

    # Интенсивность найденных линий
    lines_intensity = spectra[peak_positions_real]
    # Собираем все параметры для каждой линии в один лист
    line_nums = [num for num in range(len(NIST_wavelength))]
    parameter_list = list(zip(line_nums, NIST_wavelength, lines_intensity, NIST_A, NIST_E, NIST_G))

    parameter_list = [x for x in parameter_list if x[2] > 500]
    #parameter_list = [parameter_list[6],  parameter_list[9],  parameter_list[10]]
    peak_used_positions = [x[1] for x in parameter_list]
    peak_used_intensity = [x[2] for x in parameter_list]
    
    
    #for line in parameter_list:
        #print(line)

    k = 0.00008617333262
    temperature = []
    for num, line1 in enumerate(parameter_list):
        parameter_list.pop(num)
        for line2 in parameter_list:
            if line1[0] != line2[0] and line1[2] > 0 and line2[2] > 0:
                T = ((line1[4] - line2[4])/k) * 1/(math.log((line1[2]*line2[3]*line2[5]*line1[1]) / (line2[2]*line1[3]*line1[5]*line2[1]))) 
                if abs(T) > 0:
                    temperature.append((line1[1], line2[1], abs(T)))
                else:
                    pass
                    #print(f' false temp {(line1[1], line2[1], abs(T))}')

    return temperature
