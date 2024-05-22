# Algorithm and Architecture Optimization for 2D Discrete Fourier Transforms with Simultaneous Edge Artifact Removal

## Compilação e Execução

Para compilar o código, certifique-se de ter a biblioteca Intel MKL instalada e use o seguinte comando:

```bash
icx -o example example.c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -qopenmp