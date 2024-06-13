import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def read_cvector_bin(filename, rows, columns):
    with open(filename, "rb") as file:
        data = np.fromfile(file, dtype=np.complex64)
        if data.size != rows * columns:
            raise ValueError("Tamanho do arquivo não corresponde às dimensões da matriz")
        matrix = data.reshape((rows, columns))
    return matrix

def plot_spectrum(matrix, filename, save_path=None, cmap='viridis'):
    magnitude = np.abs(matrix)
        
    magnitude_log = np.log(magnitude + 1)
    
    plt.figure(figsize=(10, 6))
    plt.imshow(magnitude_log, cmap=cmap, aspect='auto')
    plt.colorbar(label='Magnitude (log scale)')
    plt.title(f'Magnitude Spectrum\n{filename}')
    plt.xlabel('X')
    plt.ylabel('Y')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Uso: python3 plot_complex.py <dirname> <filename> <rows> <columns> <cmap>")
        sys.exit(1)
    
    dirname = sys.argv[1]
    filename = sys.argv[2]
    rows = int(sys.argv[3])
    columns = int(sys.argv[4])
    cmap = sys.argv[5]
    filepath = f"../bin/{dirname}/{filename}.bin"
    
    if cmap not in plt.colormaps():
        print(f"Erro: '{cmap}' não é um colormap válido. Use um dos seguintes: {plt.colormaps()}")
        sys.exit(1)
    
    try:
        # Create the directory for the output image if it doesn't exist
        output_dir = f"../img/{dirname}"
        os.makedirs(output_dir, exist_ok=True)
        
        spectrum = read_cvector_bin(filepath, rows, columns)
        save_path = f"{output_dir}/{filename}.png"
        plot_spectrum(spectrum.T, filename, save_path, cmap)
    except Exception as e:
        print(f"Erro: {e}")
