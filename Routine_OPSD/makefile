# Definições
CC = gcc
CFLAGS = -fopenmp -I./include
LDFLAGS = -lmkl_rt -lm -ldl

# Diretórios
BINDIR = bin
OBJDIR = obj
SRCDIR = src
INCDIR = include

# Arquivos
EXEC = $(BINDIR)/out
SRCS = $(filter-out $(SRCDIR)/main.c, $(wildcard $(SRCDIR)/*.c))
OBJS = $(SRCS:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
MAIN = $(SRCDIR)/main.c

# Regras
all: $(EXEC)

$(EXEC): $(OBJS) $(MAIN) | $(BINDIR)
	@$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	@$(CC) $(CFLAGS) -c $< -o $@

$(BINDIR) $(OBJDIR):
	@mkdir -p $@

clean:
	@rm -rf $(OBJDIR)/*.o $(EXEC)

run: $(EXEC)
	@echo "Running with arguments: $(ARGS)"
	@$(EXEC) $(ARGS)

.PHONY: all clean
