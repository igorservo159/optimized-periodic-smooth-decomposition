import numpy as np
import matplotlib.pyplot as plt

def read_complex_matrix(filename, rows, columns):
    with open(filename, "rb") as file:
        data = np.fromfile(file, dtype=np.complex64)
        if data.size != rows * columns:
            raise ValueError("Tamanho do arquivo não corresponde às dimensões da matriz")
        matrix = data.reshape((rows, columns))
    return matrix

def plot_spectrum(matrix):
    magnitude = np.abs(matrix)
    
    # Verificar os valores da magnitude
    print(f'Magnitude Min: {magnitude.min()}, Max: {magnitude.max()}')
    
    # Aplicar transformação logarítmica
    magnitude_log = np.log(magnitude + 1)
    
    plt.imshow(magnitude_log, cmap='viridis', aspect='auto')
    plt.colorbar(label='Magnitude (log scale)')
    plt.title('Magnitude Spectrum')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

if __name__ == "__main__":
    rows, columns = 1201, 401
    spectrum = read_complex_matrix("../spectrum.bin", rows, columns)
    plot_spectrum(spectrum)
