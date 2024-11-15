# Definition
CC = gcc
CFLAGS = -fopenmp -I./include
LDFLAGS = -lmkl_rt -lm -ldl -liomp5 -lpthread

# Directories
BINDIR = bin
OBJDIR = obj
SRCDIR = src
INCDIR = include

# Files
EXEC = $(BINDIR)/out
SRCS = $(filter-out $(SRCDIR)/main.c, $(wildcard $(SRCDIR)/*.c))
OBJS = $(SRCS:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
MAIN = $(SRCDIR)/main.c

# Rules
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
