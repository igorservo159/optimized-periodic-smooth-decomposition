# Algorithm and Architecture Optimization for 2D Discrete Fourier Transforms with Simultaneous Edge Artifact Removal

Este repositório, que possui fins exclusivamente educacionais, contém a implementação de um algoritmo bastante útil no processamento digital de imagens, conhecido como Periodic Plus Smooth Decomposition (PSD). Inspirado por um artigo do International Journal of Reconfigurable Computing, esta implementação aborda a remoção de artefatos de borda em transformadas de Fourier discretas 2D, comuns em imagens não periódicas tratadas como periódicas.

O projeto fornece instruções para compilação e execução do algoritmo em C usando a biblioteca Intel MKL, além de perfilar o código com o Intel VTune para otimização de desempenho. Um exemplo de configuração JSON também é fornecido para auxiliar na configuração do ambiente de desenvolvimento.

É importante salientar que este foi inspirado por um artigo do International Journal of Reconﬁgurable Computing, da editora Hindawi, o qual está disponível no link abaixo.

[Clique aqui para acessar o artigo!](https://onlinelibrary.wiley.com/doi/10.1155/2018/1403181)

## Compilação e Execução do Algoritmo OPSD em C

Certifique-se de ter a biblioteca Intel MKL instalada!

Para compilar o código, você pode utilizar o arquivo makefile, que já está bem configurado para compilar o código, mas caso necessite, as diretivas de compilação são as seguintes:

```bash
icx -o example example.c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -qopenmp
```

Recomendo utilizar o compilador icx, da própria intel, que pode ser adquirido instalando o oneapi base toolkit, o qual já vem com várias ferramentas usadas nesse projeto, incluindo o perfilador e a MKL Intel.

## Perfilar Código com VTune

Para perfilar o código, certifique-se de ter o software instalado e use o seguinte comando:

```bash
vtune -collect hotspots -result-dir=results ./example
```

Você pode checar os resultados com:

```bash
vtune-gui results
```

Se necessário, utilize o seguinte comando para habilitar o VTune:

```bash
source /opt/intel/oneapi/vtune/latest/env/vars.sh
```

## Configuração JSON VSCODE

Abaixo está um exemplo de configuração JSON que pode ser usado no VSCode, caso você o utilize.

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

## Plotando os espectros e as imagens

Os binários das imagens e espectros de frequência podem sem encontrados no diretório bin/example/, e suas respectivas imagens no diretório img/example/.

O código C contém, em utils, funções para ler e salvar binários.

Há também, no diretório Plot/, um código em python (plot_float32.py) para plotar imagens e um para plotar os espectros de frequência (plot_complex.py). Se for utilizá-los, passe como parâmetro o nome do binário (que será também utilizado como o nome da imagem/espectro gerado), o número de rows e o número de columns como no exemplo abaixo.

```bash
python3 plot_complex.py periodic 1201 401
```