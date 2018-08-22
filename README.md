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

Now create a conda environment for `automechanic` as follows.
Make sure you have relatively recent C/C++ compilers before proceeding.
```
conda create -n amenv pip cmake python=2
source activate amenv
(amenv) conda install -c openbabel openbabel
(amenv) pip install numpy pandas future pytest
```
(The `(amenv)` indicates that you are running the command in this environment.)
The last step is only necessary if your system compilers are old.

Next, install `x2z` as follows.
```
(amenv) git clone --recurse-submodules https://github.com/PACChem/x2z  # for git version < 2.13 use --recursive
(amenv) cd x2z
(amenv) cmake . -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
(amenv) make install
```

Finally install `automechanic` into the environment.
```
(amenv) git clone https://github.com/PACChem/automechanic
(amenv) pip install automechanic/
```
You can uninstall it later with `pip uninstall automechanic`.

To test your installation, run the tests as follows.
```
(amenv) pytest automechanic/tests -v
```

## How to run the example

The example is found in `automechanic/example`.

You can parse the CHEMKIN mechanism file as follows,
where `species-in.csv` is a CSV file with columns `species` and `smiles`
listing the CHEMKIN species name and smiles string for each species in
the mechanism.
```
(amenv) automech_init -m mechanism.txt -s species-in.csv
```
This will generate geometries for all species in `geoms/` which you can
replace with your own structures if you wish.  It will also generate a file
`reactions.csv` containing SMIRKS strings for all reactions in the
mechanism, which will be classified as pressure independent, low-pressure
limit, or fall-off reactions.
The `reactions.csv` file is the input for the abstractions finder.

For testing purposes, I have provided a reduced set of sample reactions in
`reactions-sample.csv`.
To find the abstractions and initialize directories for running TorsScan, do this:
```
(amenv) automech_abstr_init -r reactions-sample.csv
```
To run abstractions 0, 1, 2, 3, 5, and 7 on node 0, do the following.
```
(amenv) automech_abstr_run -n 0 -x 0-3 5 7 -c python /path/to/torsional_scan.py
```


For a dry run, you could do something like the following instead
```
(amenv) automech_abstr_run -n 0 -x 0-3 5 7 -c cat input.dat
```
which will run `cat input` in each directory.
