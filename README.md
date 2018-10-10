# automechanic

## Installation

Clone the code:
```
git clone https://github.com/PACChem/automechanic
```

Install conda (if you haven't already):
```
export CONDA=$HOME/miniconda
export PATH=$CONDA/bin:$PATH
bash automechanic/install/conda.sh
```

Create the conda environment for `automechanic`:
```
bash automechanic/install/amenv3.sh
```
For python2, use the script for the `amenv2` environment.

Install automechanic in the environment:
```
source activate amenv3  # or amenv2
(amenv3) pip install automechanic/
```

To test your installation, run the tests as follows.
```
(amenv3) pytest automechanic/tests -v
```

To exit the environment, use `source deactivate`.


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
