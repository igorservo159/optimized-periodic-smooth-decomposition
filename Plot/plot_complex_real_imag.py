import numpy as np
import matplotlib.pyplot as plt

def read_complex_matrix(filename, rows, columns):
    with open(filename, "rb") as file:
        data = np.fromfile(file, dtype=np.complex64)
        if data.size != rows * columns:
            raise ValueError("Tamanho do arquivo não corresponde às dimensões da matriz")
        matrix = data.reshape((rows, columns))
    return matrix

def plot_matrix(matrix, title, save_path):
    plt.imshow(matrix, cmap='viridis', aspect='auto')
    plt.colorbar(label='Value')
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig(save_path)
    plt.close()

if __name__ == "__main__":
    rows, columns = 1201, 401
    spectrum = read_complex_matrix("../bin/tests/spectrum.bin", rows, columns)
    
    # Parte Real
    real_part = np.real(spectrum.T)
    plot_matrix(real_part, 'Real Part Spectrum', "../img/tests/spectrum_real.png")
    
    # Parte Imaginária
    imag_part = np.imag(spectrum.T)
    plot_matrix(imag_part, 'Imaginary Part Spectrum', "../img/tests/spectrum_imaginary.png")
