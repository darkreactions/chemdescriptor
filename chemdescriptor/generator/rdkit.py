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
                self.descriptors = whitelist['descriptors']
            else:
                raise Exception(
                    "'whitelist' should be a dict. Found: {}".format(type(whitelist)))
        else:
            self.descriptors = default_command_dict['descriptors']

        try:
            iter(input_molecules)
            self.smiles = input_molecules
        except TypeError:
            print('input_molecules is not an iterable, should be a list, tuple etc.')

        for molecule in self.smiles:
            self.molecules.append(Chem.MolFromSmiles(molecule))

    def generate(self, output_file_path, dataframe=False):
        table_data = defaultdict(list)
        table_data['Compound'] = self.smiles
        for molecule in self.molecules:
            for descriptor in self.descriptors:
                command = self.descriptors[descriptor]['command']
                column_names = self.descriptors[descriptor]['column_names']
                desc_function = self.descriptor_dict[' '.join(command)]
                table_data[' '.join(column_names)].append(
                    desc_function(molecule))
        table_df = pd.DataFrame.from_dict(table_data)
        table_df.to_csv(output_file_path, index=False)


if __name__ == "__main__":
    rd = RDKitDescriptorGenerator(
        '../examples/test_smiles.smi', '../examples/rdkit_desc_list.json')
    rd.generate('../examples/output.csv')
