# automechanic

## Installation

If you haven't already, first install some form of `conda`.
The following instructions will work for Miniconda on Linux.
```
export CONDA=$HOME/miniconda2
export PATH=$CONDA/bin:$PATH
unset PYTHONPATH
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p $CONDA
conda update conda
```

Now you can install `automechanic` as follows.
Make sure you have relatively recent C/C++ compilers before proceeding.
```
git clone https://github.com/PACChem/automechanic
git clone --recursive https://github.com/PACChem/x2z  # for git>2.13 use --recurse-submodules
conda env create -f automechanic/environment.yml
source activate amenv
(amenv) cd x2z
(amenv) cmake . -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
(amenv) make install
(amenv) cd ..
(amenv) pip install automechanic/
```

To test your installation, run the tests as follows.
```
(amenv) pytest automechanic/tests -v
```


## How to run the examples

The examples are found in `examples/` in this directory.
Here we outline the procedure for `examples/syngas`.

The `automech` program has a subcommand structure similar to programs like `git`.
You can see the available subcommands with `automech -h` or `automech --help`.

A mechanism is initialized with the `init` subcommand, which requires one of two
sets of inputs:
 1. A JSON-format mechanism file from RMG (flag: `-j`).
    Example:
    ```
    (amenv) automech init -j mechanism.json -P syngas/
    ```
 2. A CHEMKIN mechanism file (flag: `-m`).
    In this case, you must also provide a CSV file for translating your CHEMKIN
    names into meaning full species identifiers (SMILES and multiplicity).
    The CHEMKIN name should be under the column header `species` and the
    species ID should be under the column header `species_id` in the format
    `<SMILES>_m<multiplicity>`.
    Example:
    ```
    (amenv) automech init -m mechanism.txt -s species.csv -P syngas/
    ```

This initialization step will create a file called `reactions.csv` (flag: `-R`)
with all of the reactions in the mechanism.
It will also create a directory (flag: `-D`) of `.xyz` files with the
geometries for all species in the mechanism, whose paths will be written to
`species.csv` under the column header `path`.

```
(amenv) cd syngas/
(amenv) automech additions init -s species.csv -r reactions.csv -P additions/ -p
(amenv) cd additions/
# Pick reactions in reactions.csv and put them in batches/b1.csv along with a
# template input file at batches/t1.txt
(amenv) automech additions run -s ../species.csv -r reactions.csv -b batches/b1.csv -t batches/t1.txt -y nodes:d -p cmd ls
```
(The `-p` flag prints the log to screen)
