# Makefile for Matrix Multiplication Benchmark

# Detect OS
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  # On macOS M1: use Homebrew's libomp for OpenMP
  CC      := clang
  BREW_PREFIX := $(shell brew --prefix libomp)
  CFLAGS  := -g -O3 \
             -Xpreprocessor -fopenmp \
             -I$(BREW_PREFIX)/include
  IFLAGS  := -L$(BREW_PREFIX)/lib -lomp -lpthread -lm
else
  # On Linux or GCC on macOS: use AVX and OpenMP
  CC      := gcc
  CFLAGS  := -mavx -mavx2 -g -O3 -fopenmp
  IFLAGS  := -lpthread -lm
endif

SRC      := src/matmul.c src/utils.c
TARGET   := matmul

# Default parameters (can be overridden on the command line)
METHOD   ?= Sequential
BASE     ?= 64
THREADS  ?= $(shell echo $${OMP_NUM_THREADS:-1})

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(IFLAGS)

run: $(TARGET)
	@if [ -z "$(INPUT)" ]; then \
		echo "Usage: make run INPUT=<inputfile> [METHOD=<method>] [BASE=<threshold>] [THREADS=<threads>]"; \
		exit 1; \
	fi
	OMP_NUM_THREADS=$(THREADS) ./$(TARGET) -i $(INPUT) -m $(METHOD) -b $(BASE) -t $(THREADS)

experiment: $(TARGET)
	bash experiments/run_experiments.sh

clean:
	rm -f $(TARGET)

.PHONY: all run experiment clean
