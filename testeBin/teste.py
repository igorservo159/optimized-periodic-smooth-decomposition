import numpy as np
import matplotlib.pyplot as plt

# Definições das dimensões da matriz
nx, nz = 1201, 401
dx, dz = 5, 5  # Resoluções espaciais

# Ler o arquivo binário
filename = 'toy2024_cp0_5m.bin'
data = np.fromfile(filename, dtype=np.float32)

# Verificar se o tamanho dos dados está correto
if data.size != nx * nz:
    raise ValueError("Tamanho do arquivo não corresponde às dimensões da matriz")

# Redimensionar os dados para uma matriz 2D
matrix = data.reshape((nx, nz))

# Plotar o modelo
plt.figure(figsize=(10, 6))
extent = [0, (nx-1) * dx, (nz-1) * dz, 0]  # Definir os limites dos eixos
plt.imshow(matrix.T, cmap='viridis', aspect='auto', extent=extent)
plt.colorbar(label='Amplitude')
plt.title('Modelo Original')
plt.xlabel('Distância X (m)')
plt.ylabel('Profundidade Z (m)')
plt.show()
