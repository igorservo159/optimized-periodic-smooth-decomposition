import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import sys

def read_float32(filename, rows, columns):
    with open(filename, "rb") as file:
        image = np.fromfile(file, dtype=np.float32)
        if image.size != rows * columns:
            raise ValueError("Tamanho da imagem não corresponde às dimensões da matriz")
        matrix = image.reshape((rows, columns))
    return matrix

def plot_image(matrix, filename, save_path=None, cmap='grey', vmin=None, vmax=None):
    plt.figure(figsize=(10, 6))
    plt.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax,)
    plt.colorbar(label='Magnitude')
    plt.title(f'{filename}')
    plt.xlabel('X')
    plt.ylabel('Y')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def extract_center_square_window(image, window_size, rows, columns):
    rows, columns = image.shape
    center_row = rows // 2
    center_col = columns // 2
    half_n = window_size // 2
    return image[center_row - half_n:center_row + half_n, center_col - half_n:center_col + half_n]

def insert_center_square_window(image, rows, columns, window, window_size):
    rows, columns = image.shape
    center_row = rows // 2
    center_col = columns // 2
    half_n = window_size // 2

    image[center_row - half_n:center_row + half_n, center_col - half_n:center_col + half_n] = window

    return image


def apply_fft2d(image_window):
    return np.fft.fft2(image_window)

def apply_ifft2d(fft_image_window):
    return np.fft.ifft2(fft_image_window)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Uso: python3 windowing.py <dirname> <rows> <columns> <window_size> <transpose>")
        sys.exit(1)
    
    dirname = sys.argv[1]
    rows = int(sys.argv[2])
    columns = int(sys.argv[3])
    window_size = int(sys.argv[4])
    transpose = sys.argv[5].lower()
    
    dirpath = f"../bin/{dirname}"
    filepath = os.path.join(dirpath, "image.bin")

    savename = f"../bin/{dirname}/windowing"
    os.makedirs(savename, exist_ok=True)
    window_savepath = os.path.join(savename, "image.bin")
    fft_window_savepath = os.path.join(savename, "fft_window_image.bin")

    image = read_float32(filepath, rows, columns)
    if transpose == 'yes':
        image = image.T
    window = extract_center_square_window(image, window_size, rows, columns)
    window.tofile(window_savepath)
    
    #
    window_spectrum = apply_fft2d(window)
    processed_fft_window = apply_ifft2d(window_spectrum)
    fft_window = np.real(processed_fft_window)
    fft_window.astype(np.float32).tofile(fft_window_savepath)

    # Running Routine
    routine = "ccr"
    precision = "single"
    save_vector = "yes"
    input_file = "rb"
    seed = 1
    savename_ = f"{dirname}/windowing"

    args = f"{window_size} {window_size} {routine} {precision} {save_vector} {input_file} {savename_} {seed}"

    make_directory = "../Routine_OPSD"

    # Execute o comando make com os argumentos no diretório especificado
    subprocess.run(f"make run ARGS=\"{args}\"", shell=True, cwd=make_directory)

    os.rename(f"{savename}/filtered_image.bin", f"{savename}/opsd_window_image.bin")
    opsd_window_filepath = os.path.join(savename, "opsd_window_image.bin")
    
    opsd_window = read_float32(opsd_window_filepath, window_size, window_size)
    
    image_ = image.copy()
    
    reconstructed_image_by_fft = insert_center_square_window(image_, rows, columns, fft_window, window_size)
    reconstructed_image_by_opsd = insert_center_square_window(image, rows, columns, opsd_window, window_size)

    reconstructed_image_by_fft_savepath = os.path.join(dirpath, "reconstructed_image_by_fft.bin")
    reconstructed_image_by_opsd_savepath = os.path.join(dirpath, "reconstructed_image_by_opsd.bin")

    reconstructed_image_by_fft.tofile(reconstructed_image_by_fft_savepath)
    reconstructed_image_by_opsd.tofile(reconstructed_image_by_opsd_savepath)
