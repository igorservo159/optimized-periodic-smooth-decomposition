import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftshift

# Substitua 'path_to_binary_file' pelo caminho real do arquivo binário
path_to_binary_file = '../bin/tests/toy2024_cp0_5m.bin'

# Defina as dimensões da imagem (substitua pelos valores reais)
image_height = 1201  # Altura da imagem
image_width = 401   # Largura da imagem

# Leitura do arquivo binário
image_data = np.fromfile(path_to_binary_file, dtype=np.float32)
image_data = (image_data.reshape((image_height, image_width))).T

# Plot da imagem original
plt.figure(figsize=(10, 6))
plt.imshow(image_data, cmap='viridis', aspect='auto')
plt.title('Imagem Original')
plt.colorbar()
plt.show()

# Cálculo da FFT2D
fft_image = fft2(image_data)
fft_image_shifted = fftshift(fft_image)
magnitude_spectrum = np.log(np.abs(fft_image) + 1)  # Adiciona 1 para evitar log(0)

# Plot do espectro da FFT2D
plt.figure(figsize=(10, 6))
plt.imshow(magnitude_spectrum, cmap='viridis', aspect='auto')
plt.title('Espectro da FFT2D')
plt.colorbar()
plt.savefig("../img/tests/example_spectrum.png")
