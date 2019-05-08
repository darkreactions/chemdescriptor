# chemdescriptor - Molecular descriptor generator
Generic molecular descriptor generator wrapper around various software packages to simplify the process of getting descriptors

## To install
Type: ```pip install chemdescriptor```

## Requirements
1. Pandas
2. Working copy of ChemAxon cxcalc

## Usage
Currently only supports ChemAxon cxcalc. The module can be expanded to cover other generators as well

### Command Line
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

Import the generator class
``` from chemdescriptor import ChemAxonDescriptorGenerator ```

Instantiate a generator
``` cag = ChemAxonDescriptorGenerator('/path/to/SMILES/file',
                                      '/path/to/descriptor/whitelist/json',
                                      ph_values=[6, 7, 8],
                                      command_stems=None,
                                      ph_command_stems=None)
```

Generate csv output
``` cag.generate('output.csv') ```

## Notes:

Input SMILES file has a SMILES code in each line
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
These dictionaries "translate" the above descriptors into commands that chem axon can understand.

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

# To Do
[ ] Test on different machines
[ ] Get feedback on what needs to be changed/improved
[ ] Expand to cover other descriptor generators