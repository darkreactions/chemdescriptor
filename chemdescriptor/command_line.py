# from .generator import ChemAxonDescriptorGenerator, RDKitDescriptorGenerator
from .generator.chemaxon import ChemAxonDescriptorGenerator
from .generator.rdkit import RDKitDescriptorGenerator

import argparse
import json


def cxcalc():
    parser = argparse.ArgumentParser(
        prog="chemdescriptor-cx",
        description="Generate molecular descriptors from ChemAxon cxcalc",
    )
    parser.add_argument(
        "-m", "--molecule", help="Path to input SMILES file", required=True
    )
    parser.add_argument(
        "-d", "--descriptors", help="Path to descriptor white list json file"
    )
    parser.add_argument(
        "-p",
        "--pH",
        type=float,
        nargs="+",
        help="List of pH values at which to calculate descriptors",
        required=True,
    )
    parser.add_argument(
        "-c", "--commands", help="Optional command stems for descriptors in json format"
    )
    parser.add_argument("-o", "--output", help="Path to output file", required=True)
    parser.add_argument("-l", "--logfile", help="Path to a log file")

    args = parser.parse_args()
    print("Running cxcalc")
    with open(args.molecule, "r") as f:
        smiles_list = f.read().splitlines()

    whitelist = {}
    if args.descriptors:
        with open(args.descriptors, "r") as f:
            whitelist = json.load(f)

    command_dict = None
    if args.commands:
        with open(args.commands, "r") as f:
            command_dict = json.load(f)

    ph = []
    if args.pH:
        ph = args.pH

    c = ChemAxonDescriptorGenerator(
        smiles_list,
        whitelist=whitelist,
        ph_values=ph,
        command_dict=command_dict,
        logfile=args.logfile,
    )
    c.generate(args.output)


def rdkit():
    parser = argparse.ArgumentParser(
        prog="chemdescriptor-rdkit",
        description="Generate molecular descriptors from rdkit",
    )
    parser.add_argument(
        "-m", "--molecule", help="Path to input SMILES file", required=True
    )
    parser.add_argument(
        "-d", "--descriptors", help="Path to descriptor white list json file"
    )
    parser.add_argument("-o", "--output", help="Path to output file", required=True)
    parser.add_argument(
        "-c", "--commands", help="Optional command stems for descriptors in json format"
    )
    parser.add_argument("-l", "--logfile", help="Path to a log file")
    args = parser.parse_args()

    with open(args.molecule, "r") as f:
        smiles_list = f.read().splitlines()

    whitelist = {}
    if args.descriptors:
        with open(args.descriptors, "r") as f:
            whitelist = json.load(f)

    command_dict = {}
    if args.commands:
        with open(args.commands, "r") as f:
            command_dict = json.load(f)

    r = RDKitDescriptorGenerator(
        smiles_list,
        whitelist=whitelist,
        command_dict=command_dict,
        logfile=args.logfile,
    )
    r.generate(args.output)
