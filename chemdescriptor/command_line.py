from .chemaxon_generator import ChemAxonDescriptorGenerator
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--molecule',
                        help="Path to input SMILES file", required=True)
    parser.add_argument('-d', '--descriptors',
                        help="Path to descriptor white list json file", required=True)
    parser.add_argument('-p', '--pH', type=float, nargs='+',
                        help="List of pH values at which to calculate descriptors", required=True)
    parser.add_argument('-c', '--commands',
                        help="Optional command stems for descriptors in json format")
    parser.add_argument('-pc', '--phcommands',
                        help="Optional command stems for pH dependent descriptorsin json format")
    parser.add_argument('-o', '--output',
                        help='Path to output file', required=True)

    args = parser.parse_args()

    c = ChemAxonDescriptorGenerator(args.molecule,
                                    args.descriptors,
                                    ph_values=args.pH,
                                    command_stems=args.commands,
                                    ph_command_stems=args.phcommands)
    c.generate(args.output)
