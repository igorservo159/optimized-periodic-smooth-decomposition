# Algorithm and Architecture Optimization for 2D Discrete Fourier Transforms with Simultaneous Edge Artifact Removal

Este repositório, que possui fins exclusivamente educacionais, contém a implementação de um algoritmo bastante útil no processamento digital de imagens, conhecido como Periodic Plus Smooth Decomposition (PSD). Inspirado por um artigo do International Journal of Reconfigurable Computing, esta implementação aborda a remoção de artefatos de borda em transformadas de Fourier discretas 2D, comuns em imagens não periódicas tratadas como periódicas.

É importante salientar que este foi inspirado por um artigo do International Journal of Reconﬁgurable Computing, da editora Hindawi, o qual está disponível no link abaixo.

[Clique aqui para acessar o artigo!](https://onlinelibrary.wiley.com/doi/10.1155/2018/1403181)

## Compilação do Algoritmo OPSD em C

Certifique-se de ter a biblioteca Intel MKL instalada!
[Clique aqui para acessar a página de download da biblioteca](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html).

Após a instalação, lembre de utilizar as variáveis de ambiente da MKL que possuem caminho default /opt/intel/oneapi/setvars.sh.

Para compilar o código, recomendo que você utilize o arquivo makefile em Routine_OPSD, que já está bem configurado para compilar o código.
```bash
make
```

Mas caso necessite, as diretivas de compilação são as seguintes:

```bash
gcc -o example example.c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -qopenmp
```

Recomendo utilizar o compilador icx, da própria intel, que pode ser adquirido instalando o oneapi base toolkit, o qual já vem com várias ferramentas usadas nesse projeto, incluindo o perfilador e a MKL Intel. Caso queira utilizá-lo lembre de trocar o compilador de gcc para icx no makefile.

## Execução do Algoritmo OPSD em C

Quando for utilizar o algoritmo, obtenha o .bin da imagem que você quer filtrar, crie um diretório em Paper_OPSD/bin/ e coloque o binário dentro com o nome data.bin. Depois disso, é só executar o algoritmo.

Os arquivos .bin e .png são sempre salvos com os mesmos nomes ("data", "data_filtered", "spectrum", "spectrum_shifted", "smooth", "smooth_shifted", "periodic", "periodic_shifted") nos diretórios Paper_OPSD/bin/{dirname}/ e Paper_OPSD/img/{dirname}/, respectivamente. 

Para executar o algoritmo OPSD em C, é recomendado que você utilize o comando make run em Routine_OPSD forneça os seguintes argumentos: 

```bash
make run ARGS="<rows> <columns> <routine> <precision> <save_vector> <input> <dirname> <seed>"
```

Aqui estão as referências de strings para os argumentos:

- `Routine`:
  - `ccr`: compute complete routine.
  - `css`: compute shifted spectrums.
  - `cts`: compute tradicional spectrums.

- `Input`:
  - `rb`: read binary.
  - `fm`: fill matrix.
  
- `Precision`:
  - `single`
  - `double`

- `Save_vector`:
  - `yes`
  - `no`

- `Dirname`:
  - `example`
  - `other dirname in bin`


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

Os binários das imagens e espectros de frequência podem sem encontrados no diretório Paper_OPSD/bin/{dirname}/, e suas respectivas imagens no diretório Paper_OPSD/img/{dirname}/.

O código C contém, em utils, funções para ler e salvar binários.

Há também, no diretório Paper_OPSD/Plot/, um código em python (plot_float32.py) para plotar imagens e um para plotar os espectros de frequência (plot_complex.py). Se for utilizá-los, passe como parâmetro o nome do diretório do e o nome do binário (que serão também utilizados como o diretório e nome da imagem/espectro gerado), o número de rows, o número de columns, o color map e se deseja transpor a matriz como nos exemplos abaixo.

```bash
python3 plot_float.py example data 1201 401 grey yes

python3 plot_float.py example data 1201 401 viridis no

```

A imagem será salva em Paper_OPSD/img/{dirname}