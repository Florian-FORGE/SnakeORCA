## Installation instructions

## Constructing a conda env

Following the Orca installation instructions  
https://github.com/jzhoulab/orca

## Installing miniforge

If you don't have it already, follow instructions from https://github.com/conda-forge/miniforge.

### Cloning SnakeOrca

```
git clone git@forge.inrae.fr:vincent.rocher/snakeorca.git
cd snakeorca
```

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
conda deactivate
mamba env update -f snakeorca_env_part2.yml
```

Now install libstdcxx-ng==13.2.0 then pytabix :

```
conda activate snakeorca_env
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

**Ensuring access to Orca**

```
cd ..
ORCA_DIR=$(pwd)
ENV_NAME="snakeorca_env"
ENV_PATH=$(conda env list | awk -v env="$ENV_NAME" '$1 == env {print $NF}')

mkdir -p "$ENV_PATH/etc/conda/activate.d"
mkdir -p "$ENV_PATH/etc/conda/deactivate.d"

```

Add PYTHONPATH export to activation script

```
echo 'export PYTHONPATH="ORCA_DIR:$PYTHONPATH"' > "$ENV_PATH/etc/conda/activate.d/env_vars.sh"

```

Optional: unset PYTHONPATH on deactivation

```
echo 'unset PYTHONPATH' > "$ENV_PATH/etc/conda/deactivate.d/env_vars.sh"

```

**Installing snakemake**

```
conda activate snakeorca_env
mamba install -c conda-forge -c bioconda snakemake
pip install snakemake-executor-plugin-slurm

```
**Installing plotnine**

```
mamba install conda-forge::plotnine
```

### Installing SnakeOrca

```
cd ..
pip install .
```

(Ensure you are in the snakeorca directory)


## Testing the installation

```
cd orca_snake

conda activate snakeorca_env

ORCA_DIR=PATH/TO/orca

export PYTHONPATH="$ORCA_DIR":$PYTHONPATH

```

Change the the ORCA_DIR path to fit your setup.

### Testing the top-down pipeline

I you have access to genotoul, change the slurm_account in the config.yaml file in the genotoul directory. If you cannot use GPU resources change the slurm_account in the config_cpu.yaml file (rename the other one and name this one config.yaml) and change no_cuda to True in the config.yaml in the orca_snake directory.
Then use :

```
snakemake --configfile config.yaml --profile genotoul -j 3 -p -n

```

Else use :

```
snakemake --configfile config.yaml -j 3 -p -n

```

(Without GPU it can take some time to tun the test pipeline (~1h), we recommand to proceed with GPU if possible.)

**Remove the '-n' flag to run the pipeline, if the dry-run worked properly.**


### Testing the bottom-up pipeline

If the top-down pipeline properly executed, then use this to retrieve the traces from the tested mutations :


**Define the source repository and destination directory**
```
SOURCE_ROOT="work/experience/Test"
DEST_DIR="../test_resources/test_annot/bottom_up"

```

**Create destination directory if it doesn't exist**

```
mkdir -p "$DEST_DIR"

```

**Find and copy each trace_Wtd_mut.csv file**

```
find "$SOURCE_ROOT" -type f -name "trace_Wtd_mut.csv" | while read FILE; do
    REPO_NAME=$(basename "$(dirname "$(dirname "$(dirname "$FILE")")")")
    NEW_NAME="trace_${REPO_NAME}.csv"
    cp "$FILE" "$DEST_DIR/$NEW_NAME"
done

```

Then launch the second pipeline using either one of those depending on GPU accessibility :

```
snakemake --snakefile Snakefile_bottom_up --configfile config_bottom_up.yaml --profile genotoul -j 8 -p -n

```

 or :

```
snakemake --snakefile Snakefile_bottom_up --configfile config_bottom_up.yaml -j 8 -p -n

```

**Remove the '-n' flag to run the pipeline, if the dry-run worked properly.**
