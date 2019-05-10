import json
import os
import subprocess
import csv
from collections import OrderedDict
from itertools import chain
import copy
import pandas as pd

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
    version = 'x.x'
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
                 input_molecule_file_path,
                 descriptor_file_path,
                 ph_values=[],
                 command_stems=None,
                 ph_command_stems=None):
        """
        Initializes the generator
        Args:
            input_molecule_file_path:   path to ip molecules
            descriptor_file_path:       path to ip descriptors
            ph_values:                  list of pH values at which to calculate descriptors
            command_stems:              Dictonary of descriptors and its command stem
            ph_command_stems:           Dict of pH related descriptors and command stems
        """
        super().__init__(input_molecule_file_path, descriptor_file_path)

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
        with open(descriptor_file_path, 'r') as f:
            desc = json.load(f)
            self.descriptors = desc['descriptors']
            self.ph_descriptors = desc['ph_descriptors']

        # Read smiles
        with open(input_molecule_file_path, 'r') as f:
            self.smiles = f.read().splitlines()

        # Setup Descriptor commands required for chemaxon
        self._command_dict = OrderedDict()
        self._setup_descriptor_commands()
        self._setup_ph_descriptor_commands()

    def _setup_descriptor_commands(self):
        """
        Converts entries in the descriptor whitelist to the corresponding commands.
        Checks if commands are in list form, if not, convert

        """
        if self._command_stems:
            for descriptor in self.descriptors:
                if type(self._command_stems[descriptor]) is list:
                    self._command_dict[descriptor] = self._command_stems[descriptor]
                elif type(self._command_stems[descriptor]) is str:
                    self._command_dict[descriptor] = \
                        self._command_stems[descriptor].split()
                else:
                    raise Exception('Command for descriptor {} should be list or str.'
                                    'Found {}'.format(descriptor, type(self._command_stems[descriptor])))
        else:
            for cmd in self.descriptors:
                # Cannot use dict comprehension, because OrderedDict
                if type(cmd) is list:
                    self._command_dict[cmd] = cmd
                elif type(cmd) is str:
                    self._command_dict[cmd] = cmd.split()
                else:
                    raise Exception('Command for descriptor {} should be list or str.'
                                    'Found {}'.format(cmd, type(cmd)))

    def _setup_ph_descriptor_commands(self):
        """
        Generates and adds ph specific descriptors to the final
        list of descriptors.

        """

        for ph_desc in self.ph_descriptors:
            command = self._ph_command_stems[ph_desc]
            nominal_copy = copy.copy(command)
            if type(nominal_copy) is list:
                self._command_dict["{}_nominal".format(
                    ph_desc)] = nominal_copy
            elif type(nominal_copy)is str:
                self._command_dict["{}_nominal".format(ph_desc)] = \
                    nominal_copy.split()
            else:
                raise Exception('Command for descriptor {} should be list or str.'
                                'Found {}'.format(ph_desc, type(nominal_copy)))

            for pH in self.ph_values:
                ph_command = copy.copy(command)
                pH_string = str(pH).replace('.', '_')  # R compatibility
                if type(ph_command) is list:
                    ph_command.extend(['-H', str(pH)])
                elif type(ph_command) is str:
                    ph_command = ph_command.split() + ['-H', str(pH)]
                else:
                    raise Exception('Command for descriptor {} should be list or str.'
                                    'Found {}'.format(ph_desc, type(ph_command)))
                self._command_dict["{}_ph{}".format(
                    ph_desc, pH_string)] = ph_command

    def generate(self, output_file_path):
        """
        Method called to generate descriptors
        Calculates Least Energy Conformer (LEC) for each molecule by cxcalc and then writes
        into a temporary output file called lec_molecules.txt which is removed after the 
        descriptors are generated. File is written in the same folder as the desired output
        file. Inelegant solution because there is no way to directly share data with cxcalc

        Possible to use tempfile module in python?

        Args:
            output_file_path : Path to the desired output file

        """
        output_folder, output_file = os.path.split(output_file_path)
        intermediate_file = os.path.join(output_folder, 'lec_molecules.txt')
        try:
            self.generate_lec(intermediate_file)
            self.generate_descriptors(
                intermediate_file, output_file_path)
        except:
            if os.path.exists(intermediate_file):
                # Clean up intermediate file if an exception occurs
                os.remove(intermediate_file)
        os.remove(intermediate_file)

    def generate_lec(self, output_file):
        """
        Method to generate LEC. Runs cxcalc on given molecules with leconformer
        flag and outputs the LEC for each molecule

        Args:
            output_folder: Writable folder to which the intermediate file
            is written

        """

        lecProc = subprocess.run([os.path.join(self.CXCALC_PATH, 'cxcalc'), '-o',
                                  output_file,
                                  'leconformer', self.input_molecule_file_path, ])
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
            self._parse_output(output_filename)

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


if __name__ == "__main__":
    os.environ['CXCALC_PATH'] = '/home/h205c/chemaxon/bin'

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
    c = ChemAxonDescriptorGenerator('../examples/test_smile.smi',
                                    '../examples/descriptors_list.json',
                                    ph_values=[7],
                                    ph_command_stems=_cxcalcpHCommandStems)
    c.generate('./output.csv')
