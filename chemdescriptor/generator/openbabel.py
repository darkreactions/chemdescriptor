

from openbabel import pybel


molecules = list(pybel.readfile("smi", '../examples/test_foursmiles.smi'))

print(pybel.descs)
