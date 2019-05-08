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
            command_stems:              Dictonary of command stems and their name
            ph_command_stems:           Dict of pH related command stems
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
        Commands are in list form even if there is one keyword in command
        because it is easier to chain when running cxcalc

        """

        if self._command_stems:
            for descriptor in self.descriptors:
                if type(self._command_stems[descriptor]) is list:
                    self._command_dict[descriptor] = self._command_stems[descriptor]
                else:
                    self._command_dict[descriptor] = [
                        self._command_stems[descriptor]]
        else:
            for cmd in self.descriptors:
                # Cannot use dict comprehension, because OrderedDict
                self._command_dict[cmd] = [cmd]

    def _setup_ph_descriptor_commands(self):
        """
        Generates and adds ph specific descriptors to the final
        list of descriptors.
        Also generates the appropriate commands related to the pH
        descriptors
        """

        for ph_desc in self.ph_descriptors:
            command = self._ph_command_stems[ph_desc]
            if type(command) is list:
                self._command_dict["{}_nominal".format(
                    ph_desc)] = copy.copy(command)
            else:
                self._command_dict["{}_nominal".format(ph_desc)] = [
                    copy.copy(command)]

            for pH in self.ph_values:
                ph_command = copy.copy(command)
                pH_string = str(pH).replace('.', '_')  # R compatibility
                if type(ph_command) is list:
                    ph_command.extend(['-H', str(pH)])
                else:
                    ph_command = [ph_command, '-H', str(pH)]
                self._command_dict["{}_ph{}".format(
                    ph_desc, pH_string)] = ph_command

    def generate(self, output_file_path):
        """
        Method called to generate descriptors
        TODO: Find out why lec is calculated
        Args:
            output_file_path : Path to the desired output file

        """
        # output_folder, output_file = os.path.split(output_file_path)
        # lec_molecules = self.generate_lec(output_folder)
        self.generate_descriptors(
            self.input_molecule_file_path, output_file_path)

    def generate_lec(self, output_folder):
        """
        Method to filter out molecular structures with lowest
        energy conformer? Runs cxcalc on given molecules with leconformer
        flag and selects molecules whose structure is lec

        Args:
            output_folder: Writable folder to which the intermediate file
            is written

        TODO: Verify contents of output file and parse accordingly

        """
        try:
            lecProc = subprocess.run([self.CXCALC_PATH + 'cxcalc', '-o',
                                      os.path.join(
                                          output_folder, 'lec_output.txt'),
                                      'leconformer', self.input_molecule_file_path, ])
            # stdout=PIPE, stderr=PIPE, close_fds = True)
            # lec = lowest energy conformer
        except Exception:
            print("Could not run chemaxon")

    def generate_descriptors(self, smiles_molecules, output_filename):
        """
        Generate the descriptors by using cxcalc. cxcalc outputs
        to an intermediate file which is read by this script.

        Args:
            smiles_molecules:   List of SMILES molecules. Ideally
                                passed through generate_lec
            output_filename:    Path to output file

        """
        _command_list = list(chain.from_iterable(self._command_dict.values()))
        print(_command_list)
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
        data structure defined by us.
        Currently converted to pandas dataframe and written into csv.
        First column is replaced by smiles codes and column names are changed
        to whatever is defined in the input descriptors json

        Args:
            filename:       Path to output file

        """

        df = pd.read_csv(filename, sep='\t')
        df['Compound'] = self.smiles
        df = df[['Compound'] + [c for c in df if c not in ['id', 'Compound']]]

        df.columns = ['Compound'] + \
            [label for label in self._command_dict.keys()]
        df.to_csv(filename, index=False)


if __name__ == "__main__":
    os.environ['CXCALC_PATH'] = '/home/h205c/chemaxon/bin'
    """
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
    c = ChemAxonDescriptorGenerator('test_smile.smi',
                                    'descriptors_list.json',
                                    ph_values=[7],
                                    ph_command_stems=_cxcalcpHCommandStems)
    c.generate('./output.csv')
    """
