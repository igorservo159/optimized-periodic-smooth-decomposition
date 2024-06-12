import numpy as np
import matplotlib.pyplot as plt
import sys

def read_cvector_bin(filename, rows, columns):
    with open(filename, "rb") as file:
        data = np.fromfile(file, dtype=np.complex64)
        if data.size != rows * columns:
            raise ValueError("Tamanho do arquivo não corresponde às dimensões da matriz")
        matrix = data.reshape((rows, columns))
    return matrix

def plot_spectrum(matrix, filename, save_path=None):
    magnitude = np.abs(matrix)
        
    magnitude_log = np.log(magnitude + 1)
    
    plt.figure(figsize=(10, 6))
    plt.imshow(magnitude_log, cmap='viridis', aspect='auto')
    plt.colorbar(label='Magnitude (log scale)')
    plt.title(f'Magnitude Spectrum\n{filename}')
    plt.xlabel('X')
    plt.ylabel('Y')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python3 plot_complex.py <dirname> <filename> <rows> <columns>")
        sys.exit(1)
    
    dirname = sys.argv[1]
    filename = sys.argv[2]
    rows = int(sys.argv[3])
    columns = int(sys.argv[4])
    filepath = f"../bin/{dirname}/{filename}.bin"
    
    try:
        spectrum = read_cvector_bin(filepath, rows, columns)
        save_path = f"../img/{dirname}/{filename}.png"
        plot_spectrum(spectrum.T, filename, save_path)
    except Exception as e:
        print(f"Erro: {e}")
