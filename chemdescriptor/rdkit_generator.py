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

from .generator import BaseDescriptorGenerator


class RDKitDescriptorGenerator(BaseDescriptorGenerator):
    name = 'RDKit Descriptor Generator'
    __version__ = '0.0.1'

    def __init__(self,
                 input_molecules,
                 descriptors,
                 ph_values=[],
                 command_stems=None,
                 ph_command_stems=None):
        """Handles only SMILES and stores into self.molecules
        TODO: Handle InchI
        """
        super().__init__()

        self.descriptor_dict = {desc[0]: desc[1]
                                for desc in Descriptors.descList}
        self.descriptor_list_3d = [desc for desc in dir(
            Descriptors3D) if '__' not in desc]

        self.molecules = []

        # print(self.descriptor_dict.keys())

        if isinstance(descriptors, str):
            with open(descriptors, 'r') as f:
                desc = json.load(f)
        elif isinstance(descriptors, dict):
            desc = descriptors
        else:
            raise Exception(
                "'descriptors' should be a path or dict. Found: {}".format(type(descriptors)))

        self.descriptors = desc['descriptors']

        if isinstance(input_molecules, str):
            if os.path.exists(input_molecules):
                with open(input_molecules, 'r') as f:
                    self.smiles = f.read().splitlines()
        elif isinstance(input_molecules, (list, set, tuple)):
            self.smiles = input_molecules
        else:
            raise Exception(
                "'input_molecules' should be a path or list. Found: {}".format(type(input_molecules)))

        for molecule in self.smiles:
            self.molecules.append(Chem.MolFromSmiles(molecule))

    def generate(self, output_file_path, dataframe=False, lec=False):
        table_data = defaultdict(list)
        table_data['Compound'] = self.smiles
        for molecule in self.molecules:
            for descriptor in self.descriptors:
                # print(descriptor)
                desc_function = self.descriptor_dict[descriptor]
                table_data[descriptor].append(desc_function(molecule))
        table_df = pd.DataFrame.from_dict(table_data)
        table_df.to_csv(output_file_path, index=False)


if __name__ == "__main__":
    rd = RDKitDescriptorGenerator(
        '../examples/test_smiles.smi', '../examples/rdkit_desc_list.json')
    rd.generate('../examples/output.csv')
