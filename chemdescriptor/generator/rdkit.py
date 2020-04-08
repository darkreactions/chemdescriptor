import json
import os
import subprocess
import csv
from collections import OrderedDict, defaultdict
from itertools import chain
import copy
import pandas as pd
import traceback
from rdkit import Chem
from rdkit.Chem import Descriptors, Descriptors3D

from .base import BaseDescriptorGenerator
from ..defaults.rdkit import default_command_dict


class RDKitDescriptorGenerator(BaseDescriptorGenerator):
    name = 'RDKit Descriptor Generator'
    __version__ = '0.0.2'

    def __init__(self,
                 input_molecules,
                 whitelist={},
                 command_dict={},
                 logfile=None):
        """Handles only SMILES and stores into self.molecules
        TODO: Handle InchI
        """
        super().__init__(logfile=logfile)

        self.descriptor_dict = {desc[0]: desc[1]
                                for desc in Descriptors.descList}
        self.descriptor_list_3d = [desc for desc in dir(
            Descriptors3D) if '__' not in desc]

        self.molecules = []

        # print(self.descriptor_dict.keys())

        if whitelist:
            if isinstance(whitelist, dict):
                self.descriptor_whitelist = whitelist['descriptors']
            else:
                raise Exception(
                    "'whitelist' should be a dict. Found: {}".format(type(whitelist)))
        else:
            self.descriptor_whitelist = list(
                default_command_dict['descriptors'].keys())

        if command_dict:
            if isinstance(command_dict, dict):
                self.command_dict = command_dict
            else:
                raise Exception(
                    "'command_dict' should be a dict. Found: {}".format(type(command_dict)))
        else:
            self.command_dict = default_command_dict['descriptors']

        try:
            iter(input_molecules)
            self.smiles = input_molecules
        except TypeError:
            print('input_molecules is not an iterable, should be a list, tuple etc.')

        # for molecule in self.smiles:
        #    self.molecules.append(Chem.MolFromSmiles(molecule))

    def generate(self, output_file_path, dataframe=False):
        table_data = defaultdict(list)

        # for molecule in self.molecules:
        for smi in self.smiles:
            molecule = Chem.MolFromSmiles(smi)
            if molecule:
                table_data['Compound'].append(smi)
                for descriptor in self.descriptor_whitelist:
                    if descriptor not in self.command_dict:
                        self.logger.error(
                            'Descriptor {} not found in command dict'.format(descriptor))
                        continue
                    else:
                        command = self.command_dict[descriptor]['command']
                        column_names = self.command_dict[descriptor]['column_names']
                        # TODO: Join with space or underscore?
                        desc_function = self.descriptor_dict[' '.join(command)]
                        table_data[' '.join(column_names)].append(
                            desc_function(molecule))
            else:
                self.logger.error(
                    '{} did not generate a mol object'.format(smi))
        table_df = pd.DataFrame.from_dict(table_data)

        if output_file_path:
            table_df.to_csv(output_file_path, index=False)

        return table_df


if __name__ == "__main__":
    rd = RDKitDescriptorGenerator(
        '../examples/test_smiles.smi', '../examples/rdkit_desc_list.json')
    rd.generate('../examples/output.csv')
