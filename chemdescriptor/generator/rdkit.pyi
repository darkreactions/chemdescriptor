from typing_extensions import Protocol
import pandas as pd

class RDKitDescriptorGenerator(Protocol):
    def __init__(
        self,
        input_molecules: "list[str]",
        whitelist: "dict[str, str]",
        command_dict: "dict[str, str]",
        logfile: str,
    ) -> None: ...
    def generate(
        self, output_file_path: str, dataframe: bool
    ) -> None | pd.DataFrame: ...
