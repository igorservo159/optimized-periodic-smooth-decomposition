# Algorithm and Architecture Optimization for 2D Discrete Fourier Transforms with Simultaneous Edge Artifact Removal

## Compilação e Execução

Para compilar o código, certifique-se de ter a biblioteca Intel MKL instalada e use o seguinte comando:

```bash
icx -o example example.c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -qopenmp
```

## Perfilar Código com VTune

Para perfilar o código, certifique-se de ter o software instalado e use o seguinte comando:

```bash
vtune -collect hotspots -result-dir=results ./example
```

Você pode checar os resultados com:

```bash
vtune-gui results
```

Se precisar, utilize o seguinte comando para habilitar o VTune:

```bash
source /opt/intel/oneapi/vtune/latest/env/vars.sh
```


## Configuração JSON

Abaixo está um exemplo de configuração JSON que pode ser usada no seu ambiente.

```json
{
    "configurations": [
        {
            "name": "Linux",
            "includePath": [
                "${default}",
                "${workspaceFolder}/**",
                "/opt/intel/oneapi/compiler/latest/include"
            ],
            "defines": [],
            "compilerPath": "/opt/intel/oneapi/compiler/latest/bin/icx",
            "cStandard": "c17",
            "cppStandard": "c++17",
            "intelliSenseMode": "linux-clang-x64"
        }
    ],
    "version": 4
}
```

[Clique aqui para acessar o artigo que explica mais sobre o algoritmo e que embasa este projeto!](https://onlinelibrary.wiley.com/doi/10.1155/2018/1403181)
