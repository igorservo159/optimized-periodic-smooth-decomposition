# Variáveis
CC = icx
CFLAGS = -qopenmp -I./include
LDFLAGS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# Nome do executável
TARGET = main

# Diretório dos arquivos fonte
SRC_DIR = src

# Arquivos fonte
SRCS = $(SRC_DIR)/main.c $(SRC_DIR)/fourier.c $(SRC_DIR)/utils.c

# Alvo padrão
all: $(TARGET)

# Regra para compilar o alvo
$(TARGET): $(SRCS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

# Limpeza dos arquivos gerados
clean:
	rm -f $(TARGET)
