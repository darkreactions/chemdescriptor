# chemdescriptor - Molecular descriptor generator
Generic molecular descriptor generator wrapper around various software packages to simplify the process of getting descriptors

## To install
Type:

```pip install chemdescriptor```

## Requirements
1. Pandas
2. Working copy of ChemAxon cxcalc

## Usage
Currently only supports ChemAxon cxcalc. The module can be expanded to cover other generators as well.
Example input files can be found in the examples/ folder of this repo as well as the pip installed package.

**Important! The code requires an environment variable CXCALC_PATH to be set, which points to the folder where cxcalc is installed!**

### Command Line
```
chemdescriptor-cx -m /path/to/SMILES/file -d /path/to/descriptor/whitelist/json -p 6.8 7.0 7.2 -o output.csv
```

```
usage: chemdescriptor-cx [-h] -m MOLECULE -d DESCRIPTORS -p PH [PH ...]
                         [-c COMMANDS] [-pc PHCOMMANDS] -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -m MOLECULE, --molecule MOLECULE
                        Path to input SMILES file
  -d DESCRIPTORS, --descriptors DESCRIPTORS
                        Path to descriptor white list json file
  -p PH [PH ...], --pH PH [PH ...]
                        List of pH values at which to calculate descriptors
  -c COMMANDS, --commands COMMANDS
                        Optional command stems for descriptors in json format
  -pc PHCOMMANDS, --phcommands PHCOMMANDS
                        Optional command stems for pH dependent descriptorsin
                        json format
  -o OUTPUT, --output OUTPUT
                        Path to output file
```

### In code
Set CXCALC_PATH

```
import os
os.environ['CXCALC_PATH'] = '/path/to/cxcalc'
```

Import the generator class

``` from chemdescriptor import ChemAxonDescriptorGenerator ```

Instantiate a generator
``` 
cag = ChemAxonDescriptorGenerator('/path/to/SMILES/file',
                                  '/path/to/descriptor/whitelist/json',
                                  ph_values=[6, 7, 8],
                                  command_stems=None,
                                  ph_command_stems=None)
```

Generate csv output
``` cag.generate('output.csv', dataframe=False, lec=False) ```

Optional keyword arguments for `generate` include `dataframe` boolean (default False) which returns a pandas dataframe in addition to writing a csv if True
and `lec` boolean (default False) which converts the Smiles code to an intermediate "Low Energy Conformer (LEC)" representation before generating descriptors.
A license is most likely required to generate LECs.

## Notes:

Input SMILES file has a SMILES code in each line.

Descriptor whitelist is a json file of the form:
```
{
    "descriptors": [
        "refractivity",
        "maximalprojectionarea",
        "maximalprojectionradius",
        "maximalprojectionsize",
        "minimalprojectionarea",
        "minimalprojectionradius",
        "minimalprojectionsize"
    ],
    "ph_descriptors": [
        "avgpol",
        "molpol",
        "vanderwaals",
        "asa",
        "asa+",
        "asa-",
        "asa_hydrophobic",
        "asa_polar",
        "hbda_acc",
        "hbda_don",
        "polar_surface_area"
    ]
}
```

chemdescriptor expects 2 keys where "descriptors" are generic and "ph_descriptors" are ph dependent descriptors

2 optional dictionaries can be passed to the ChemAxonDescriptorGenerator, "command_stems" and "ph_command_stems".
These dictionaries "translate" the above descriptors into commands that ChemAxon cxcalc can understand.

For example, if no value is passed to the ph_command_stems, the following dictionary is used:

```
_default_ph_command_stems = {
        'avgpol': 'avgpol',
        'molpol': 'molpol',
        'vanderwaals': 'vdwsa',
        'asa': ['molecularsurfacearea', '-t', 'ASA'],
        'asa+': ['molecularsurfacearea', '-t', 'ASA+'],
        'asa-': ['molecularsurfacearea', '-t', 'ASA-'],
        'asa_hydrophobic': ['molecularsurfacearea', '-t', 'ASA_H'],
        'asa_polar': ['molecularsurfacearea', '-t', 'ASA_P'],
        'hbda_acc': 'acceptorcount',
        'hbda_don': 'donorcount',
        'polar_surface_area': 'polarsurfacearea',
    }
```

Note that commands with multiple words are entries in a list. For example, the command 

```molecularsurfacearea -t ASA```

is represented in the dictionary as ```['molecularsurfacearea', '-t', 'ASA']```

# To Do
[ ] Test on different machines

[ ] Get feedback on what needs to be changed/improved

[ ] Expand to cover other descriptor generators
