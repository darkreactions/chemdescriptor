from abc import ABCMeta, abstractmethod


import json


class BaseDescriptorGenerator(metaclass=ABCMeta):
    name = 'Base Generator'
    version = '0.0'

    def __init__(self, input_molecule_file_path, descriptor_file_path):
        self.input_molecule_file_path = input_molecule_file_path
        self.descriptor_file_path = descriptor_file_path

    @abstractmethod
    def generate(self):
        pass
