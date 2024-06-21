import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def read_float32(filename, rows, columns):
    with open(filename, "rb") as file:
        data = np.fromfile(file, dtype=np.float32)
        if data.size != rows * columns:
            raise ValueError("Tamanho do arquivo não corresponde às dimensões da matriz")
        matrix = data.reshape((rows, columns))
    return matrix

def plot_image(matrix, filename, save_path=None, cmap='viridis'):    
    plt.figure(figsize=(10, 6))
    plt.imshow(matrix, cmap=cmap, aspect='auto')
    plt.colorbar(label='Amplitude')
    plt.title(f'Imagem de {filename}.bin')
    plt.xlabel('Distância X (m)')
    plt.ylabel('Profundidade Z (m)')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Uso: python3 plot_data.py <dirname> <rows> <columns> <cmap> <transpose>")
        sys.exit(1)
    
    dirname = sys.argv[1]
    filename = "data"
    rows = int(sys.argv[2])
    columns = int(sys.argv[3])
    cmap = sys.argv[4]
    transpose = sys.argv[5].lower()
    filepath = f"../bin/{dirname}/{filename}.bin"
    
    if cmap not in plt.colormaps():
        print(f"Erro: '{cmap}' não é um colormap válido. Use um dos seguintes: {plt.colormaps()}")
        sys.exit(1)
    
    try:
        # Create the directory for the output image if it doesn't exist
        output_dir = f"../img/{dirname}"
        os.makedirs(output_dir, exist_ok=True)

        spectrum = read_float32(filepath, rows, columns)
        
        if transpose == 'yes':
            spectrum = spectrum.T
        
        save_path = f"{output_dir}/{filename}.png"
        plot_image(spectrum, filename, save_path, cmap)
    except Exception as e:
        print(f"Erro: {e}")
