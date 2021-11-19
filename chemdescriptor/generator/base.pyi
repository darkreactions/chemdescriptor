from typing_extensions import Protocol

class BaseDescriptorGenerator(Protocol):
    def __init__(self, logfile: str) -> None: ...
