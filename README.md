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
To run a given set of abstractions with TorsScan you can use the `automech_abstr_run`
command.
A dry run would look as follows
```
(amenv) automech_abstr_run -n b444 b445 -x 0-3 5 7 cmd ls -la
```
which will create inputs for abstractions 0, 1, 2, 3, 5, and 7 to run on nodes b444 and b445,
enter each job directory and run the command `ls -la`.
Anything that comes after `cmd` will be run as a command in the job directory, so an actual
run would look as follows.
```
(amenv) automech_abstr_run -n <node list> -x <abstraction id ranges> cmd python /path/to/torsional_scan.py
```
