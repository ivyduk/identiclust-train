# identiclust-train
this repo contains a workflow built using snakeMake Tool and its goal is to train a model with whole genome sequences to classify between two different species    

# requirements
* Python ≥3.3
* Snakemake 3.11.0

The easiest way to setup these prerequisites is to use the Miniconda Python 3 distribution. The tutorial assumes that you are using either Linux or MacOS X. Both Snakemake and Miniconda work also under Windows, but the Windows shell is too different to be able to provide generic examples.

## Install miniconda on Linux
``` 
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ bash Miniconda3-latest-Linux-x86_64.sh
```

## Install miniconda on MacOSX
```
$ curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
$ bash Miniconda3-latest-MacOSX-x86_64.sh
```

miniconda 3

## Create virtual environment for project and activate

```
$ conda env create --name identiclust-env --file environment.yaml
$ source activate identiclust-env
```

if your require a new package you can install using 


```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda install blast
```

## Run the workflow
```
$ snakemake
```
