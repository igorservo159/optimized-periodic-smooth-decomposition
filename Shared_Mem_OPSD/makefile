# Variáveis
CC = icx
CFLAGS = -qopenmp
LDFLAGS = -lmkl_rt -lm -ldl

# Nome do executável
TARGET = main

# Arquivos fonte
SRCS = main.c aux.c

# Alvo padrão
all: $(TARGET)

# Regra para compilar o alvo
$(TARGET): $(SRCS)
	$(CC) -o $@ $^ $(LDFLAGS) $(CFLAGS)

# Limpeza dos arquivos gerados
clean:
	rm -f $(TARGET)
