import numpy as np
import matplotlib.pyplot as plt

def read_float32(filename, rows, columns):
    with open(filename, "rb") as file:
        data = np.fromfile(file, dtype=np.float32)
        if data.size != rows * columns:
            raise ValueError("Tamanho da imagem não corresponde às dimensões da matriz")
        matrix = data.reshape((rows, columns))
    return matrix

def plot_spectrum(matrix, filename, save_path=None):
    magnitude = np.abs(matrix)
    magnitude_log = np.log(magnitude + 1)
    
    plt.figure(figsize=(10, 6))
    plt.imshow(magnitude_log.T, cmap='viridis', aspect='auto')
    plt.colorbar(label='Magnitude (log scale)')
    plt.title(f'Magnitude Spectrum')
    plt.xlabel('Frequência X')
    plt.ylabel('Frequência Y')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def fft2d_spectrum(filename, rows, columns, save_path=None):
    # Read the float32 binary file
    matrix = read_float32(filename, rows, columns)
    
    # Compute the 2D FFT and shift the zero frequency component to the center
    fft_matrix = np.fft.fft2(matrix)
    fft_matrix = np.fft.fftshift(fft_matrix)
    
    # Plot the magnitude spectrum
    plot_spectrum(fft_matrix, filename, save_path)

# Example usage
if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) != 4:
        print("Uso: python fft2d_spectrum.py <dirname> <rows> <columns>")
        sys.exit(1)
    
    dirname = sys.argv[1]
    filename = "data"
    rows = int(sys.argv[2])
    columns = int(sys.argv[3])
    filename_ = f"{filename}.bin"
    dirname_ = f"../bin/{dirname}/"

    try:
        # Determine the save path for the output image
        output_dir = "../img/testes/"
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(dirname))[0]}.png")
        readpath = os.path.join(dirname_, filename_)
        # Calculate and plot the 2D FFT spectrum
        fft2d_spectrum(readpath, rows, columns, save_path)
        print(f"Salvou espectro de frequência: {save_path}")
    except Exception as e:
        print(f"Erro: {e}")
