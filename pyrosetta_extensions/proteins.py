from pathlib import Path
from typing import Union
from collections import OrderedDict

from Bio.PDB import PDBList
from functional import seq
from pyrosetta import *
from pyrosetta.rosetta.core.pose import Pose, PDBInfo


class Protein:
    '''
    Protein class that makes working with proteins easier
    '''
    def __init__(self, path: Union[str, Path, Pose], folder_to_download = "/tmp"):
        if type(path) == Pose:
            self.pose: Pose = path
        else:
            path_str = path if type(path) == str else path.as_posix()
            if Path(path_str).exists():
                self.pose: Pose = pose_from_pdb(path_str)
            else:
                print(f"cannot find file {path_str}, trying to download it to{folder_to_download}")
                pdbl = PDBList()
                native_pdb = pdbl.retrieve_pdb_file(path_str, pdir=folder_to_download, file_format='pdb')
                self.pose: Pose = pose_from_pdb(native_pdb)

    @property
    def residues(self) -> seq:
        return seq(self.pose.residues)

    @property
    def protein_residues(self) -> seq:
        return self.residues.filter(lambda r:  r.is_protein())

    def residue(self, num: int) -> pyrosetta.rosetta.core.conformation.Residue:
        return self.pose.residue(num)

    @property
    def first(self) -> pyrosetta.rosetta.core.conformation.Residue:
        return self.residues[0]

    @property
    def info(self) -> PDBInfo:
        return self.pose.pdb_info()

    @property
    def conformation(self) -> pyrosetta.rosetta.core.conformation.Conformation:
        return self.pose.conformation()

    def pdb_range(self, start: int, end: int) -> seq:
        return seq([self.info.pose2pdb(i) for i in range(start, end)])

    def with_deleted_residues(self, start: int, end: int):
        novel_pose = self.pose.clone()
        novel_pose.delete_residue_range_slow(start, end)
        return Protein(novel_pose)

    @property
    def annotated_sequence(self) -> str:
        return self.pose.annotated_sequence()


    @property
    def sequence(self) -> str:
        return self.pose.sequence()

    @property
    def by_chain(self) -> OrderedDict:
        from collections import OrderedDict
        proteins = seq(self.pose.clone().split_by_chain()).map(lambda p: Protein(p))
        dic = OrderedDict()
        for p in proteins:
            letter = p.info.pose2pdb(1).split(" ")[1]
            dic[letter] = p
        return dic

    @property
    def by_chain_sequences(self) -> OrderedDict:
        proteins = seq(self.pose.clone().split_by_chain()).map(lambda p: Protein(p))
        dic = OrderedDict()
        for p in proteins:
            letter = p.info.pose2pdb(1).split(" ")[1]
            dic[letter] = p.sequence
        return dic
