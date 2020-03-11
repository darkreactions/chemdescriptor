# chemdescriptor - Molecular descriptor generator
Generic molecular descriptor generator wrapper around various software packages to simplify the process of getting descriptors

## To install
Type:

```pip install chemdescriptor```

## Requirements
1. Pandas
2. ChemAxon descriptors
    - Working copy of ChemAxon cxcalc
3. RDKit descriptors
    - RDKit installed

## Usage
Currently only supports ChemAxon cxcalc and RDKit. The module can be expanded to cover other generators as well.
Example input files can be found in the examples/ folder of this repo as well as the pip installed package.

### CXCalc

**Important! The code requires an environment variable CXCALC_PATH to be set, which points to the folder where cxcalc is installed!**

### Command Line
```
chemdescriptor-cx -m /path/to/SMILES/file -d /path/to/descriptor/whitelist/json -p 6.8 7.0 7.2 -o output.csv
```

```
usage: chemdescriptor-cx [-h] -m MOLECULE -d DESCRIPTORS -p PH [PH ...]
                         [-c COMMANDS] -o OUTPUT

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
  -o OUTPUT, --output OUTPUT
                        Path to output file
```

### In code

The package will initially search cxcalc executable in the PATH variable if not
will fall back to CXCALC_PATH

Set CXCALC_PATH

```
import os
os.environ['CXCALC_PATH'] = '/path/to/cxcalc'
```

Import the generator class

``` from chemdescriptor.generator.chemaxon import ChemAxonDescriptorGenerator as CAG```

Import SMILES and whitelist

```
with open('/path/to/SMILES/file', 'r') as f:
    smiles_list = f.read().splitlines()

with open('/path/to/descriptor/whitelist/json', 'r') as f:
    whitelist = json.load(f)
```

Instantiate a generator. ```smiles_list``` is a list of smiles and ```whitelist```
is a dictionary of keys in the command_dict 
```logfile``` is the path to a log which contains information such as the final cxcalc
command, columns that were renamed and other errors for debugging

``` 
cag = CAG(smiles_list,
          whitelist,
          ph_values=[6, 7, 8],
          command_dict={},
          logfile='/path/to/logfile')
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

An optional dictionary can be passed to the ChemAxonDescriptorGenerator, "command_dict" which
"translates" the above descriptors into commands that ChemAxon cxcalc can understand.

It also consists of column names that will be added to the final output

An example of a command_dict is:

```
command_dict = {
    "descriptors": {
        "atomcount_c": {
            "command": [
                "atomcount",
                "-z",
                "6"
            ],
            "column_names": [
                "_feat_AtomCount_C"
            ]
        },
        "wateraccessiblesurfacearea": {
            "command": [
                "wateraccessiblesurfacearea"
            ],
            "column_names": [
                "_feat_ASA",
                "_feat_ASA+",
                "_feat_ASA-",
                "_feat_ASA_H",
                "_feat_ASA_P"
            ]
        }
    "ph_descriptors": {
        "acceptorcount": {
            "command": [
                "acceptorcount"
            ],
            "column_names": [
                "_feat_Hacceptorcount"
            ]
        },
        "donorcount": {
            "command": [
                "donorcount"
            ],
            "column_names": [
                "_feat_Hdonorcount"
            ]
        }
    }

```
```command_dict``` consists of 2 dictionaries with keys ```descriptors``` and 
```ph_descriptors```. Within each dictionary are descriptors referred in the whitelist. 

Under each descriptor, two lists are required ```command``` and ```column_names```

Command refers to the command line options for cxcalc as documented 
[here](https://docs.chemaxon.com/display/docs/cxcalc+calculator+functions)
Note that commands with multiple words are entries in a list. For example, the command 
```atomcount -z 6``` is represented in the dictionary as ```['atomcount', '-z', '6']```

```column_names``` is a list of names the user wants to rename the cxcalc generated
csv column names.

Certain commands generate multiple columns for example, ```wateraccessiblesurfacearea```
generates 5 columns. Therefore, the ```column_names``` list becomes
```
"column_names": [
                "_feat_ASA",
                "_feat_ASA+",
                "_feat_ASA-",
                "_feat_ASA_H",
                "_feat_ASA_P"
            ]
```

*Note* : If the number of columns generated by cxcalc do not match the expected count, 
none of the column names are renamed.

### RDKit

Much easier to use. Only needs a list of descriptors similar to cxcalc. 


# To Do
[ ] Test on different machines

[ ] Get feedback on what needs to be changed/improved

[ ] Expand to cover other descriptor generators
