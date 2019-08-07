import json
import os
import subprocess
import csv
from collections import OrderedDict
from itertools import chain
import copy
import pandas as pd
import traceback

from .generator import BaseDescriptorGenerator


class ChemAxonDescriptorGenerator(BaseDescriptorGenerator):
    """
    This descriptor generator class requires a valid installation
    of ChemAxon.

    TODO: Specify output type (csv, database etc.)
    TODO: Specify input type (currently expects a random
          assortment of csv and text files)

    """
    name = 'ChemAxon Descriptor Generator'
    version = '0.0.3'
    CXCALC_PATH = None
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

    def __init__(self,
                 input_molecules,
                 descriptors,
                 ph_values=[],
                 command_stems=None,
                 ph_command_stems=None):
        """
        Initializes the generator
        Args:
            input_molecules:   path to ip molecules
            descriptors:       path to ip descriptors
            ph_values:                  list of pH values at which to calculate descriptors
            command_stems:              Dictonary of descriptors and its command stem
            ph_command_stems:           Dict of pH related descriptors and command stems
        """
        super().__init__(input_molecules, descriptors)

        # Look for cxcalc path
        if 'CXCALC_PATH' in os.environ:
            self.CXCALC_PATH = os.environ['CXCALC_PATH']
        else:
            raise Exception('CXCALC_PATH environment variable not set')

        # Initialize inputs
        self.ph_values = ph_values
        self._command_stems = command_stems
        if ph_command_stems:
            self._ph_command_stems = ph_command_stems
        else:
            self._ph_command_stems = self._default_ph_command_stems

        # Read descriptors form given json file
        if isinstance(descriptors, str):
            with open(descriptors, 'r') as f:
                desc = json.load(f)
        elif isinstance(descriptors, dict):
            desc = descriptors
        else:
            raise Exception(
                "'descriptors' should be a path or dict. Found: {}".format(type(descriptors)))

        self.descriptors = desc['descriptors']
        self.ph_descriptors = desc['ph_descriptors']

        # Read smiles
        if isinstance(input_molecules, str):
            with open(input_molecules, 'r') as f:
                self.smiles = f.read().splitlines()
        elif isinstance(input_molecules, (list, set, tuple)):
            self.smiles = input_molecules

        else:
            raise Exception(
                "'input_molecules' should be a path or list. Found: {}".format(type(input_molecules)))

        # Setup Descriptor commands required for chemaxon
        self._command_dict = OrderedDict()
        self._setup_descriptor_commands()
        self._setup_ph_descriptor_commands()

    def _setup_descriptor_commands(self):
        """
        Converts entries in the descriptor whitelist to the corresponding commands.
        Checks if commands are in list form, if not, convert

        """
        for descriptor in self.descriptors:
            self._command_dict = self._add_command(descriptor, self._command_stems,
                                                   self._command_dict)

    def _setup_ph_descriptor_commands(self):
        """
        Generates and adds ph specific descriptors to the final
        list of descriptors.

        """

        for ph_desc in self.ph_descriptors:
            self._command_dict = self._add_command(ph_desc, self._ph_command_stems,
                                                   self._command_dict, postfix='nominal')

            for ph in self.ph_values:
                ph_string = str(ph).replace('.', '_')  # R compatibility
                self._command_dict = self._add_command(ph_desc, self._ph_command_stems,
                                                       self._command_dict, postfix='ph'+ph_string, ph=ph)

    def _add_command(self, descriptor, command_stems, command_dict, postfix=None, ph=None):
        """
        Adds a specific command to the command_dict. Checks whether descriptor exists in
        the command_stems. Adds user defined postfix to command_dict key and ph related commands

        Args:
            descriptor:     Name of the descriptor being added
            command_stems:  Lookup dict to get commands for a descriptor
            command_dict:   Aggregation of all commands required for this run
            postfix:        User defined extention of the key in command_dict (Optional)
            ph:             pH value at which to run command (Optional)

        Return:
            command_dict:   Dict of all commands + current command
        """
        if command_stems and descriptor in command_stems:
            cmd = command_stems[descriptor]
        else:
            cmd = descriptor

        cmd = copy.copy(cmd)

        if postfix:
            descriptor_key = "{}_{}".format(descriptor, postfix)
        else:
            descriptor_key = descriptor

        if isinstance(cmd, list):
            command_dict[descriptor_key] = cmd
        elif isinstance(cmd, str):
            command_dict[descriptor_key] = cmd.split()
        else:
            raise Exception('Command for descriptor {} should be list or str.'
                            'Found {}'.format(descriptor, type(cmd)))

        if ph:
            command_dict[descriptor_key] += ['-H', str(ph)]

        return command_dict

    def generate(self, output_file_path, dataframe=False, lec_conformer=False):
        """
        Method called to generate descriptors
        Calculates Least Energy Conformer (LEC) for each molecule by cxcalc and then writes
        into a temporary output file called lec_molecules.txt which is removed after the 
        descriptors are generated. File is written in the same folder as the desired output
        file. Inelegant solution because there is no way to directly share data with cxcalc

        Possible to use tempfile module in python?

        Args:
            output_file_path : Path to the desired output file
            dataframe : Bool value to return a dataframe or csv file
            lec_conformer : Bool to indicate whether lec_conformer is available, default is False

        """
        output_folder, output_file = os.path.split(output_file_path)

        try:
            if lec_conformer:
                intermediate_file = os.path.join(
                    output_folder, 'lec_molecules.txt')
                self.generate_lec(intermediate_file)
            else:
                intermediate_file = self.input_molecule_file_path

            result_dataframe = self.generate_descriptors(
                intermediate_file, output_file_path)
        except Exception as e:
            print("Exception : {}".format(e))

            if lec_conformer and os.path.exists(intermediate_file):
                # Clean up intermediate file if an exception occurs
                os.remove(intermediate_file)

        if lec_conformer and os.path.exists(intermediate_file):
            os.remove(intermediate_file)

        if dataframe:
            return result_dataframe
        else:
            return None

    def generate_lec(self, output_file):
        """
        Method to generate LEC. Runs cxcalc on given molecules with leconformer
        flag and outputs the LEC for each molecule

        Args:
            output_folder: Writable folder to which the intermediate file
            is written

        """
        folder, filename = os.path.split(output_file)
        input_molecule_file_path = os.path.join(folder, 'input_smiles.smi')
        # print(self.smiles)
        with open(input_molecule_file_path, 'w') as f:
            f.writelines("\n".join(self.smiles))
        lecProc = subprocess.run([os.path.join(self.CXCALC_PATH, 'cxcalc'), '-o',
                                  output_file,
                                  'leconformer', input_molecule_file_path, ])
        os.remove(input_molecule_file_path)
        if lecProc.returncode != 0:
            print(lecProc.stderr)

    def generate_descriptors(self, smiles_molecules, output_filename):
        """
        Generate the descriptors by using cxcalc. cxcalc outputs
        to an intermediate file which is read by this script.

        Args:
            smiles_molecules:   Path to file of SMILES molecules or LEC molecules
                                Ideally passed through generate_lec
            output_filename:    Path to output file

        """
        _command_list = list(chain.from_iterable(self._command_dict.values()))

        calcProc = subprocess.run([os.path.join(self.CXCALC_PATH, 'cxcalc'),
                                   '-o', output_filename, smiles_molecules] +
                                  _command_list)
        if calcProc.returncode != 0:
            print(calcProc.stderr)
        else:
            return self._parse_output(output_filename)

    def _parse_output(self, filename):
        """
        Parses the output file generated by cxcalc and stored in a
        data structure defined by user.
        Currently converted to pandas dataframe and written into csv.
        First column is replaced by smiles codes and column names are changed
        to whatever is defined in the input descriptors json

        Args:
            filename:       Path to output file

        """

        df = pd.read_csv(filename, sep='\t')
        df['Compound'] = self.smiles
        # Replace First column with Compound
        df = df[['Compound'] + [c for c in df if c not in ['id', 'Compound']]]
        # Replace cxcalc generated labels with the ones defined in descriptor dict
        df.columns = ['Compound'] + \
            [label for label in self._command_dict.keys()]
        df.to_csv(filename, index=False)

        return df


if __name__ == "__main__":
    os.environ['CXCALC_PATH'] = '/Applications/MarvinSuite/bin'

    _cxcalcpHCommandStems = {
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
    c = ChemAxonDescriptorGenerator('../examples/test_smiles.smi',
                                    '../examples/descriptors_list.json',
                                    ph_values=[7],
                                    ph_command_stems=_cxcalcpHCommandStems)
    c.generate('output.csv')
