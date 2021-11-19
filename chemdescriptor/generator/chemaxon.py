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
from pathlib import Path

from .base import BaseDescriptorGenerator
from ..defaults.cxcalc import default_command_dict


class ChemAxonDescriptorGenerator(BaseDescriptorGenerator):
    """
    This descriptor generator class requires a valid installation
    of ChemAxon.
    """

    name = "ChemAxon Descriptor Generator"
    __version__ = "0.0.6"

    def __init__(
        self,
        input_molecules,
        whitelist={},
        ph_values=[],
        command_dict=None,
        logfile=None,
        standardize=False,
    ):
        """
        Initializes the generator
        Args:
            input_molecules:    list of smiles
            whitelist:          dict of keys in the command dict. If empty, use all descriptors
            ph_values:          list of pH values at which to calculate descriptors
            command_dict:       Dictionary of commands to run. Uses internal dict if none provided
            logfile:            Path to logfile in case of errors
        """
        super().__init__(logfile)
        self.logger.info("Starting cxcalc")
        # Look for cxcalc path
        if shutil.which("cxcalc"):
            self.CXCALC_PATH = Path(shutil.which("cxcalc")).parent
        elif "CXCALC_PATH" in os.environ:
            self.CXCALC_PATH = os.environ["CXCALC_PATH"]
        else:
            raise Exception(
                "cxcalc command not found or CXCALC_PATH environment variable not set"
            )

        # Initialize inputs

        self.ph_values = ph_values

        if command_dict:
            self._command_dict = command_dict
        else:
            self._command_dict = default_command_dict

        self.descriptors_dict = self._command_dict["descriptors"]
        self.ph_descriptors_dict = self._command_dict["ph_descriptors"]

        try:
            iter(input_molecules)
            if standardize:
                if shutil.which("standardize"):
                    self.STANDARDIZE_PATH = Path(shutil.which("standardize")).parent
                elif "STANDARDIZE_PATH" in os.environ:
                    self.STANDARDIZE_PATH = os.environ["STANDARDIZE_PATH"]
                else:
                    raise Exception(
                        "standardize command not found or STANDARDIZE_PATH environment variable not set"
                    )
                self.smiles = self._standardize(input_molecules)
            else:
                self.smiles = input_molecules
        except TypeError:
            print("input_molecules is not an iterable, should be a list, tuple etc.")

        # If the whitelist is empty, use all default descriptors
        self.whitelist = {}
        if whitelist:
            self.whitelist = whitelist
        else:
            self.whitelist["descriptors"] = list(self.descriptors_dict.keys())
            self.whitelist["ph_descriptors"] = list(self.ph_descriptors_dict.keys())

        # Setup Descriptor commands required for chemaxon
        self._command_dict = OrderedDict()
        self._column_names = []
        self._setup_descriptor_commands()
        self._setup_ph_descriptor_commands()

    def _standardize(self, smiles):
        stdProc = subprocess.run(
            [
                os.path.join(self.STANDARDIZE_PATH, "standardize"),
                *smiles,
                "-c",
                "removefragment:method=keeplargest",
            ],
            stdout=subprocess.PIPE,
        )
        smiles = [line.decode("utf-8") for line in stdProc.stdout.splitlines()]
        self.logger.info("List of standardized smiles: {}".format(smiles))
        return smiles

    def _setup_descriptor_commands(self):
        """
        Converts entries in the descriptor whitelist to the corresponding commands.
        Checks if commands are in list form, if not, convert

        """
        for descriptor in self.whitelist["descriptors"]:
            self._add_command(descriptor, self.descriptors_dict)

    def _setup_ph_descriptor_commands(self):
        """
        Generates and adds ph specific descriptors to the final
        list of descriptors.

        """

        for ph_desc in self.whitelist["ph_descriptors"]:
            self._add_command(ph_desc, self.ph_descriptors_dict, postfix="", ph=3.0)

            for ph in self.ph_values:
                ph_string = str(ph).replace(".", "_")  # R compatibility
                self._add_command(
                    ph_desc, self.ph_descriptors_dict, postfix="ph" + ph_string, ph=ph
                )

    def _add_command(self, descriptor, command_stems, postfix=None, ph=None):
        """
        Adds a specific command to the command_dict. Checks whether descriptor exists in
        the command_stems. Adds user defined postfix to command_dict key and ph related commands

        Args:
            descriptor:     Name of the descriptor being added
            command_stems:  Lookup dict to get commands for a descriptor
            command_dict:   Aggregation of all commands required for this run
            postfix:        User defined extention of the key in command_dict (Optional)
            ph:             pH value at which to run command (Optional)

        """

        if command_stems and (descriptor in command_stems):
            cmd = command_stems[descriptor]["command"]
            col_names = command_stems[descriptor]["column_names"]
        else:
            cmd = descriptor
            col_names = descriptor

        cmd = copy.copy(cmd)
        col_names = copy.copy(col_names)

        if postfix:
            descriptor_key = "{}_{}".format(descriptor, postfix)
            col_names = ["{}_{}".format(name, postfix) for name in col_names]
        else:
            descriptor_key = descriptor

        if isinstance(cmd, list):
            self._command_dict[descriptor_key] = cmd
            self._column_names += col_names
        else:
            raise Exception(
                "Command for descriptor {} should be list"
                "Found {}".format(descriptor, type(cmd))
            )

        if ph:
            self._command_dict[descriptor_key] += ["-H", str(ph)]

    def generate(
        self, output_file_path, dataframe=False, lec=False, rename_columns=False
    ):
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

        self.input_molecule_file_path = os.path.join(output_folder, "input_smiles.smi")

        with open(self.input_molecule_file_path, "w") as f:
            f.writelines("\n".join(self.smiles))
        try:
            if lec:
                intermediate_file = os.path.join(output_folder, "lec_molecules.txt")
                self.generate_lec(intermediate_file)
                result_dataframe = self.generate_descriptors(
                    intermediate_file, output_file_path, rename_columns
                )
            else:
                result_dataframe = self.generate_descriptors(
                    self.input_molecule_file_path, output_file_path, rename_columns
                )
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

        lecProc = subprocess.run(
            [
                os.path.join(self.CXCALC_PATH, "cxcalc"),
                "-o",
                output_file,
                "leconformer",
                self.input_molecule_file_path,
            ]
        )
        if lecProc.returncode != 0:
            self.logger.error(lecProc.stderr)
            print(lecProc.stderr)

    def generate_descriptors(self, smiles_molecules, output_filename, rename_columns):
        """
        Generate the descriptors by using cxcalc. cxcalc outputs
        to an intermediate file which is read by this script.

        Args:
            smiles_molecules:   Path to file of SMILES molecules or LEC molecules
                                Ideally passed through generate_lec OR a list of smiles strings
            output_filename:    Path to output file

        """
        _command_list = list(chain.from_iterable(self._command_dict.values()))
        self.logger.info("Command list: {}".format(_command_list))

        calcProc = subprocess.run(
            [
                os.path.join(self.CXCALC_PATH, "cxcalc"),
                "-g",
                "-o",
                output_filename,
                smiles_molecules,
            ]
            + _command_list
        )
        if calcProc.returncode != 0:
            print(calcProc.stderr)
            self.logger.error(calcProc.stderr)

        return self._parse_output(output_filename, rename_columns)

    def _parse_output(self, filename, rename_columns):
        """
        Parses the output file generated by cxcalc and stored in a
        data structure defined by user.
        Currently converted to pandas dataframe and written into csv.
        First column is replaced by smiles codes and column names are changed
        to whatever is defined in the input descriptors json

        Args:
            filename:       Path to output file

        """

        df = pd.read_csv(filename, sep="\t")
        df["Compound"] = self.smiles
        # Replace First column with Compound
        df = df[["Compound"] + [c for c in df if c not in ["id", "Compound"]]]

        if rename_columns:
            # Replace cxcalc generated labels with the ones defined in descriptor dict
            # Commented out because cxcalc does not generate the same number of columns
            if len(self._column_names) + 1 == len(df.columns):
                for row in zip(df.columns[1:], self._column_names):
                    self.logger.info(row)
                df.columns = ["Compound"] + [label for label in self._column_names]

            else:
                self.logger.error("Final column names: {}".format(self._column_names))
                self.logger.error("Current column names: {}".format(df.columns))

        df.to_csv(filename, index=False)

        return df


if __name__ == "__main__":
    os.environ["CXCALC_PATH"] = "C:/Program Files/ChemAxon/MarvinSuite/bin"
    os.environ[
        "STANDARDIZE_PATH"
    ] = "C:/Program Files/ChemAxon/StructureRepresentationToolkit/bin"

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

    with open("../../examples/cxcalc_command_dict.json", "r") as f:
        command_dict = json.load(f)

    with open("../../examples/test_smiles_short.smi", "r") as f:
        smiles = f.read().splitlines()

    whitelist = {
        "descriptors": [key for key in default_command_dict["descriptors"].keys()],
        "ph_descriptors": [
            key for key in default_command_dict["ph_descriptors"].keys()
        ],
    }

    cag = ChemAxonDescriptorGenerator(
        smiles,
        whitelist={},
        ph_values=[6.1, 7],
        command_dict=command_dict,
        logfile="output.log",
        standardize=True,
    )
    cag.generate("output.csv", lec=False)
