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


## Instructions for plotting:

To parse Arrhenius parameters from a CHEMKIN file:
```
automech chemkin to_csv mechanism.txt -p
```
For our purposes, let's add the flags `-R ck_reactions.csv -S ck_species.csv`
to rename the outputs.
To get the reaction IDs you must have a `species.csv` file with the species IDs
for each CHEMKIN species name, in which case you can run
```
automech chemkin id_reactions ck_reactions.csv species.csv -R ck_reactions.csv -p
```
to add a `reaction_id` column to `ck_reactions.csv`.

To extract Arrhenius parameters from the output files after a run:
```
automech reactions find_arrhenius reactions.csv -p
```
where `reactions.csv` contains the paths to the run directories.
You may then want to perform a sort on one of the Arrhenius coefficient
columns:
```
automech csv sort reactions.csv arrh_a -p
```
which will move the non-empty rows to the top of the file.

To write the Arrhenius parameters in CHEMKIN format:
```
automech reactions to_chemkin reactions.csv -p
```
which will generate a file `mechanism.txt`.

To plot the Arrhenius parameters:
```
automech reactions plot_arrhenius reactions.csv -p
```
which will create a directory of image files and write their paths to
`reactions.csv`.
To open the files, you may need to install an image viewer:
```
pip install pygame image-view
```
after which you can open the images using `image-view /path/to/file`.

More elaborate Arrhenius plots can also be generated:
```
automech reactions plot_arrhenius reactions.csv -r path/to/ck_reactions.csv -k reaction reaction_id -x 100 2000 -e png -p
```
This will plot the calculated rate coefficient (blue) against the one found in
the CHEMKIN file (orange), which we wrote to `ck_reactions.csv` above, and print
both the reaction ID and the CHEMKIN reaction name on the figure.
