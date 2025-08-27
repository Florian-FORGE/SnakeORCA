## Installation instructions

## Constructing a conda env

Following the Orca installation instructions  
https://github.com/jzhoulab/orca

### Installing the conda snake orca env

```
mamba env create -f snakeorca_env_part1.yml
mamba activate snakeorca_env
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
```

Log on a gpu node and make sure cuda is available
```
srun --mem 16G -p gpuq --gres=gpu:A100_1g.10gb:1 --pty bash
python -c 'import torch; print(torch.cuda.is_available())'
```
The last command should return True

```
mamba deactivate
mamba env update -f snakeorca_env_part2.yml
```

Now install libstdcxx-ng==13.2.0 then pytabix :

```
mamba activate snakeorca_env
mamba install conda-forge::libstdcxx-ng==13.2.0
mamba install bioconda::pytabix
```

**Installing Orca**

```
git clone https://github.com/jzhoulab/orca.git
cd orca
wget https://zenodo.org/record/6234936/files/resources_core.tar.gz
wget https://zenodo.org/record/6234936/files/resources_mcools.tar.gz
wget https://zenodo.org/record/4594676/files/resources_extra.tar.gz
tar xf resources_core.tar.gz
tar xf resources_mcools.tar.gz
tar xf resources_extra.tar.gz

git clone https://github.com/kathyxchen/selene.git
cd selene
git checkout custom_target_support
python setup.py build_ext --inplace
python setup.py install
```

**Installing snakemake**

```
mamba activate snakeorca_env
mamba install -c conda-forge -c bioconda snakemake
pip install snakemake-executor-plugin-slurm

```
**Installing plotnine**

```
mamba install conda-forge::plotnine
```

### Installing SnakeOrca

```
git clone git@forge.inrae.fr:vincent.rocher/snakeorca.git
cd snakeorca
pip install .
```



## Testing the installation

```
cd orca_snake

conda activate snakeorca_env

ORCA_DIR=PATH/TO/orca

export PYTHONPATH="$ORCA_DIR":$PYTHONPATH

snakemake --configfile config.yaml --profile genotoul -j 3 -p -n

```
