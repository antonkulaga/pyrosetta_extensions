import click
from openmm.app import *
from openmm import *
from openmm.unit import *
from pyrosetta import *
from pyrosetta.rosetta import *
from pathlib import Path
from pyrosetta.rosetta.core.pose import Pose



@click.command()
@click.option('--pdb', help='pdb load', default="/data/for_docking/positive_controls/1adq_fixed.pdb")
def load(pdb: str):
    print(f"pdb is {pdb}")
    pose: Pose = pose_from_pdb(pdb)
    print(f"total residues: {pose.total_residue()}")
    for i in pose.residues:
        print(i)

if __name__ == '__main__':
    pyrosetta.init()
    load()