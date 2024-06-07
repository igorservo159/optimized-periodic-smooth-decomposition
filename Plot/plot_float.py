import numpy as np
import matplotlib.pyplot as plt
import sys

def read_float32(filename, rows, columns):
    with open(filename, "rb") as file:
        data = np.fromfile(file, dtype=np.float32)
        if data.size != rows * columns:
            raise ValueError("Tamanho do arquivo não corresponde às dimensões da matriz")
        matrix = data.reshape((rows, columns))
    return matrix

def plot_image(matrix, filename, save_path=None):    
    plt.figure(figsize=(10, 6))
    plt.imshow(matrix.T, cmap='viridis', aspect='auto')
    plt.colorbar(label='Amplitude')
    plt.title(f'Imagem de {filename}.bin')
    plt.xlabel('Distância X (m)')
    plt.ylabel('Profundidade Z (m')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python3 plot_float.py <nome_do_binario> <rows> <columns>")
        sys.exit(1)
    
    filename = sys.argv[1]
    rows = int(sys.argv[2])
    columns = int(sys.argv[3])
    filepath = f"../bin/example/{filename}.bin"
    
    try:
        spectrum = read_float32(filepath, rows, columns)
        save_path = f"../img/example/{filename}.png"
        plot_image(spectrum, filename, save_path)
    except Exception as e:
        print(f"Erro: {e}")
