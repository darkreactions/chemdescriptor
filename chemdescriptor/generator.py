from abc import ABCMeta, abstractmethod


import json


class BaseDescriptorGenerator(metaclass=ABCMeta):
    name = 'Base Generator'
    version = '0.0'

    def __init__(self):
        pass

    @abstractmethod
    def generate(self):
        pass
