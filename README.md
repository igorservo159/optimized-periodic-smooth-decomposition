# Algorithm and Architecture Optimization for 2D Discrete Fourier Transforms with Simultaneous Edge Artifact Removal

Este repositório, que possui fins exclusivamente educacionais, contém a implementação de um algoritmo bastante útil no processamento digital de imagens, conhecido como Periodic Plus Smooth Decomposition (PSD). Inspirado por um artigo do International Journal of Reconfigurable Computing, esta implementação aborda a remoção de artefatos de borda em transformadas de Fourier discretas 2D, comuns em imagens não periódicas tratadas como periódicas.

É importante salientar que este foi inspirado por um artigo do International Journal of Reconﬁgurable Computing, da editora Hindawi, o qual está disponível no link abaixo.

[Clique aqui para acessar o artigo!](https://onlinelibrary.wiley.com/doi/10.1155/2018/1403181)

## Compilação do Algoritmo OPSD em C

Certifique-se de ter a biblioteca Intel MKL instalada!

Para compilar o código, recomendo que você utilize o arquivo makefile, que já está bem configurado para compilar o código, mas caso necessite, as diretivas de compilação são as seguintes:

```bash
icx -o example example.c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -qopenmp
```

Recomendo utilizar o compilador icx, da própria intel, que pode ser adquirido instalando o oneapi base toolkit, o qual já vem com várias ferramentas usadas nesse projeto, incluindo o perfilador e a MKL Intel.

## Execução do Algoritmo OPSD em C

Para executar o algoritmo OPSD em C, utilize o executável bin/out e forneça os seguintes argumentos: `<rows> <columns> <routine> <precision> <save_vector> <input> <seed>`.

Aqui estão as referências de strings para os argumentos:

- `Routine`:
  - `ccr`: rotina de cálculo concluída.
  - `css`: espectros deslocados computados.
  - `cts`: espectros tradicionais computados.

- `Precision`:
  - `single`: precisão de ponto flutuante de precisão simples.
  - `double`: precisão de ponto flutuante dupla.

- `Save_vector`:
  - `yes`: salvar vetor.
  - `no`: não salvar vetor.

- `Input`:
  - `rb`: ler arquivo binário.
  - `fm`: preencher matriz.

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

## Plotando os espectros e as imagens

Os binários das imagens e espectros de frequência podem sem encontrados no diretório bin/example/, e suas respectivas imagens no diretório img/example/.

O código C contém, em utils, funções para ler e salvar binários.

Há também, no diretório Plot/, um código em python (plot_float32.py) para plotar imagens e um para plotar os espectros de frequência (plot_complex.py). Se for utilizá-los, passe como parâmetro o nome do binário (que será também utilizado como o nome da imagem/espectro gerado), o número de rows e o número de columns como no exemplo abaixo.

```bash
python3 plot_complex.py periodic 1201 401
```