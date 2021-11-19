from typing_extensions import Protocol
import pandas as pd

class ChemAxonDescriptorGenerator(Protocol):
    def __init__(
        self,
        input_molecules: "list[str]",
        whitelist: "dict[str, str]",
        ph_values: "list[float]",
        command_dict: "dict[str, str]",
        logfile: str,
    ) -> None: ...
    def generate(
        self, output_file_path: str, dataframe: bool, lec: bool
    ) -> None | pd.DataFrame: ...
    def generate_lec(self, output_file: str): ...
    def generate_descriptors(
        self, smiles_molecules: str | "list[str]", output_filename: str
    ): ...
