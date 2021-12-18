#Python
from typing import Union
from Bio.PDB import PDBList
from pyrosetta import *
from pyrosetta.rosetta import *
from pathlib import Path
from pyrosetta.rosetta.core.pose import Pose, PDBInfo
from functional import seq
from pyrosetta.rosetta.protocols import antibody


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
    def first(self) -> pyrosetta.rosetta.core.conformation.Residue:
        return self.residues[0]

    @property
    def info(self) -> PDBInfo:
        return self.pose.pdb_info()

    @property
    def conformation(self) -> pyrosetta.rosetta.core.conformation.Conformation:
        return self.pose.conformation()


class Antibody(Protein):

    def __init__(self, path: Union[str, Path, Pose], numbering = antibody.Chothia, folder_to_download="/tmp"):
        Protein.__init__(self, path, folder_to_download)
        self.antibody_info: antibody.AntibodyInfo = antibody.AntibodyInfo(self.pose, antibody.AHO_Scheme, numbering)

    @property
    def antibody_sequence(self) -> str:
        return self.antibody_info.get_antibody_sequence()

    @property
    def numbering_schema(self) -> str:
        return  self.antibody_info.get_current_AntibodyNumberingScheme()

    @property
    def something(self):
        self.antibody_info.get_AntibodyFrameworkInfo()

    @property
    def antibody_chain_names(self) -> str:
        return self.antibody_info.get_antibody_chain_string()

    @property
    def start_H1(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.H1, self.pose)

    @property
    def start_H2(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.H2, self.pose)

    @property
    def start_H3(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.H3, self.pose)

    @property
    def start_H4(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.H4, self.pose)

    @property
    def sequence_H1(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.H1, self.pose)

    @property
    def sequence_H2(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.H2, self.pose)

    @property
    def sequence_H3(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.H3, self.pose)


    @property
    def sequence_H4(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.H4, self.pose)


    @property
    def start_L1(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.L1, self.pose)

    @property
    def start_L2(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.L2, self.pose)

    @property
    def start_L3(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.L3, self.pose)

    @property
    def start_L4(self):
        return self.antibody_info.get_CDR_start(antibody.CDRNameEnum.L4, self.pose)

    @property
    def sequence_L1(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.L1, self.pose)

    @property
    def sequence_L2(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.L2, self.pose)

    @property
    def sequence_L3(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.L3, self.pose)

    @property
    def sequence_L4(self):
        return self.antibody_info.get_CDR_sequence_with_stem(antibody.CDRNameEnum.L4, self.pose)

