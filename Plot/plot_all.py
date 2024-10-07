import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def read_float32(filename, rows, columns):
    with open(filename, "rb") as file:
        image = np.fromfile(file, dtype=np.float32)
        if image.size != rows * columns:
            raise ValueError("Tamanho da imagem não corresponde às dimensões da matriz")
        matrix = image.reshape((rows, columns))
    return matrix

def read_cvector_bin(filename, rows, columns):
    with open(filename, "rb") as file:
        spectrum = np.fromfile(file, dtype=np.complex64)
        if spectrum.size != rows * columns:
            raise ValueError("Tamanho do espectro não corresponde às dimensões da matriz")
        matrix = spectrum.reshape((rows, columns))
    return matrix

def plot_image(matrix, filename, save_path=None, cmap='viridis', vmin=None, vmax=None):
    plt.figure(figsize=(10, 6))
    plt.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax,)
    #plt.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, extent=[0, matrix.shape[1]*5, matrix.shape[0]*5, 0])
    plt.colorbar(label='Magnitude')
    plt.title(f'{filename}')
    plt.xlabel('X')
    plt.ylabel('Y')
    #plt.xlabel('Distância X (m)')
    #plt.ylabel('Profundidade Z (m)')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_spectrum(matrix, filename, save_path=None):
    magnitude = np.abs(matrix)
    magnitude_log = np.log(magnitude + 1)
    
    plt.figure(figsize=(10, 6))
    plt.imshow(magnitude_log, cmap='viridis', aspect='auto')
    plt.colorbar(label='Magnitude (log scale)')
    plt.title(f'{filename}')
    plt.xlabel('X')
    plt.ylabel('Y')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Uso: python3 plot_all.py <dirname> <rows> <columns> <cmap> <transpose>")
        sys.exit(1)
    
    dirname = sys.argv[1]
    rows = int(sys.argv[2])
    columns = int(sys.argv[3])
    cmap = sys.argv[4]
    transpose = sys.argv[5].lower()
    dirpath = f"../bin/{dirname}"
    
    if cmap not in plt.colormaps():
        print(f"Erro: '{cmap}' não é um colormap válido. Use um dos seguintes: {plt.colormaps()}")
        sys.exit(1)
    
    try:
        # Create the directory for the output images if it doesn't exist
        output_dir = f"../img/{dirname}"
        os.makedirs(output_dir, exist_ok=True)
        
        # Read the image.bin file to determine vmin and vmax
        image_filepath = os.path.join(dirpath, "image.bin")
        image_matrix = read_float32(image_filepath, rows, columns)
        vmin, vmax = image_matrix.min(), image_matrix.max()
        
        # List all binary files in the directory
        for filename in os.listdir(dirpath):
            if filename.endswith(".bin"):
                filepath = os.path.join(dirpath, filename)
                
                try:
                    # Determine the type of data and read accordingly
                    if any(x in filename for x in ["spectrum"]):
                        spectrum = read_cvector_bin(filepath, rows, columns)
                        if transpose == 'yes':
                            spectrum = spectrum.T
                        save_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.png")
                        plot_spectrum(spectrum, filename, save_path)
                    elif any (x in filename for x in ["image"]):
                        image = read_float32(filepath, rows, columns)
                        if transpose == 'yes':
                            image = image.T
                        save_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.png")
                        plot_image(image, filename, save_path, cmap, vmin, vmax)
                    
                    print(f"Salvou imagem: {save_path}")
                except Exception as e:
                    print(f"Erro ao processar {filename}: {e}")
    except Exception as e:
        print(f"Erro: {e}")
