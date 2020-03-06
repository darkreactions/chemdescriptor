import json
import os
import subprocess
import csv
from collections import OrderedDict
from itertools import chain
import copy
import pandas as pd
import traceback
import shutil

from chemdescriptor.generator import BaseDescriptorGenerator
from chemdescriptor.cxcalc_defaults import default_command_dict


class ChemAxonDescriptorGenerator(BaseDescriptorGenerator):
    """
    This descriptor generator class requires a valid installation
    of ChemAxon.

    TODO: Specify output type (csv, database etc.)
    TODO: Specify input type (currently expects a random
          assortment of csv and text files)

    """
    name = 'ChemAxon Descriptor Generator'
    __version__ = '0.0.6'

    def __init__(self,
                 input_molecules,
                 whitelist=[],
                 ph_values=[],
                 command_dict=None,
                 logfile=None):
        """
        Initializes the generator
        Args:
            input_molecules:    list of smiles
            whitelist:          list of keys in the command dict. If empty, use all descriptors
            ph_values:          list of pH values at which to calculate descriptors
            command_dict:       Dictionary of commands to run. Uses internal dict if none provided
            logfile:            Path to logfile in case of errors
        """
        super().__init__(logfile)

        # Look for cxcalc path
        # TODO: Auto search for cxcalc
        print(shutil.which('cxcalc'))
        if 'CXCALC_PATH' in os.environ:
            self.CXCALC_PATH = os.environ['CXCALC_PATH']
        else:
            raise Exception('CXCALC_PATH environment variable not set')

        # Initialize inputs
        self.ph_values = ph_values

        if command_dict:
            self._command_dict = command_dict
        else:
            self._command_dict = default_command_dict

        self.descriptors = self._command_dict['descriptors']
        self.ph_descriptors = self._command_dict['ph_descriptors']

        try:
            iter(input_molecules)
            self.smiles = input_molecules
        except TypeError:
            print('input_molecules is not an iterable, should be a list, tuple etc.')

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
        #print('descriptor: {}'.format(descriptor))
        #print('command_stems: {}'.format(command_stems))
        if command_stems and (descriptor in command_stems):
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

    def generate(self, output_file_path, dataframe=False, lec=False):
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
            lec : Bool to indicate whether leconformer is available in cxcalc, default is False

        """
        output_folder, output_file = os.path.split(output_file_path)

        self.input_molecule_file_path = os.path.join(output_folder,
                                                     'input_smiles.smi')

        with open(self.input_molecule_file_path, 'w') as f:
            f.writelines("\n".join(self.smiles))
        try:
            if lec:
                intermediate_file = os.path.join(
                    output_folder, 'lec_molecules.txt')
                self.generate_lec(intermediate_file)
                result_dataframe = self.generate_descriptors(
                    intermediate_file, output_file_path)
            else:
                result_dataframe = self.generate_descriptors(
                    self.input_molecule_file_path, output_file_path)
        except Exception as e:
            print("Exception : {}".format(e))

            if lec and os.path.exists(intermediate_file):
                # Clean up intermediate file if an exception occurs
                os.remove(intermediate_file)

        if lec and os.path.exists(intermediate_file):
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
                                Ideally passed through generate_lec OR a list of smiles strings
            output_filename:    Path to output file

        """
        _command_list = list(chain.from_iterable(self._command_dict.values()))
        print(_command_list)

        calcProc = subprocess.run([os.path.join(self.CXCALC_PATH, 'cxcalc'),
                                   '-g', '-o', output_filename, smiles_molecules] +
                                  _command_list)
        if calcProc.returncode != 0:
            # print(calcProc.stderr)
            self.logger.error(calcProc.stderr)

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
        # Commented out because cxcalc does not generate the same number of columns
        if len(self._command_dict.keys()) + 1 == len(df.columns):
            df.columns = ['Compound'] + \
                [label for label in self._command_dict.keys()]
        else:
            print(self._command_dict.keys())
            print(df.columns)

        df.to_csv(filename, index=False)

        return df


if __name__ == "__main__":
    os.environ['CXCALC_PATH'] = '/Applications/MarvinSuite/bin'

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
    c = ChemAxonDescriptorGenerator('/Users/vshekar/Code/misc_test_code/test_smiles.txt',
                                    '/Users/vshekar/Code/misc_test_code/descriptors_list.json',
                                    ph_values=[])
    c.generate('output.csv', lec=True)
    
    command_stems = json.load(open(
        '/Users/vshekar/Downloads/debug_for_shekar/descriptor_files/descriptor_commands_set1_v0.json', 'r'))
    ph_command_stems = json.load(open(
        '/Users/vshekar/Downloads/debug_for_shekar/descriptor_files/descriptor_ph_commands_set1_v0.json', 'r'))
    """

    with open('../examples/cxcalc_command_dict.json', 'r') as f:
        command_dict = json.load(f)

    with open('../examples/test_smiles.smi', 'r') as f:
        smiles = f.read().splitlines()

    whitelist = []

    cag = ChemAxonDescriptorGenerator(smiles,
                                      whitelist=whitelist,
                                      ph_values=[6],
                                      command_dict=command_dict)
    cag.generate('output.csv', lec=False)
