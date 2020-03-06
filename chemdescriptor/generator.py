from abc import ABCMeta, abstractmethod
import logging

import json


class BaseDescriptorGenerator(metaclass=ABCMeta):
    name = 'Base Generator'
    version = '0.0'

    def __init__(self, logfile=None):
        self.logger = logging.getLogger('chemdescriptor')
        if logfile:
            hdlr = logging.FileHandler(logfile)
            formatter = logging.Formatter(
                '%(asctime)s %(levelname)s %(message)s')
            hdlr.setFormatter(formatter)
            logger.addHandler(hdlr)

    @abstractmethod
    def generate(self):
        pass
