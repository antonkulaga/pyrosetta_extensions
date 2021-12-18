from pathlib import Path
from typing import Union
from collections import OrderedDict

from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.protocols import antibody

from pyrosetta_extensions.proteins import Protein





class Antibody(Protein):

    def __init__(self, path: Union[str, Path, Pose], numbering = antibody.Chothia, folder_to_download="/tmp"):
        Protein.__init__(self, path, folder_to_download)
        self.antibody_info: antibody.AntibodyInfo = antibody.AntibodyInfo(self.pose, antibody.AHO_Scheme, numbering)
        self.enum_manager = antibody.AntibodyEnumManager()

    @property
    def antibody_sequence(self) -> str:
        self.antibody_info.qq_heavy_residue()
        return self.antibody_info.get_antibody_sequence()

    @property
    def numbering_schema(self) -> str:
        return  self.antibody_info.get_current_AntibodyNumberingScheme()

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

    def with_deleted_residues(self, start: int, end: int) -> 'Antibody':
        novel_pose = self.pose.clone()
        novel_pose.delete_residue_range_slow(start, end)
        return Antibody(novel_pose)


    def print_regions(self):
        for i in range(1, self.pose.size()+1):
            region = self.antibody_info.get_region_of_residue(self.pose, i)
            residue = self.residue(i)
            if region == antibody.cdr_region:
                region_str = self.enum_manager.cdr_name_enum_to_string(self.antibody_info.get_CDRNameEnum_of_residue(self.pose, i))
            else:
                region_str = self.enum_manager.antibody_region_enum_to_string(region)
            print(f"{i} {self.info.pose2pdb(i)} {residue.name1()} {residue.name()} {region_str}")


    @property
    def end_H1(self):
        return self.antibody_info.get_CDR_end(antibody.CDRNameEnum.H1, self.pose)

    @property
    def end_H2(self):
        return self.antibody_info.get_CDR_end(antibody.CDRNameEnum.H2, self.pose)

    @property
    def end_H3(self):
        return self.antibody_info.get_CDR_end(antibody.CDRNameEnum.H3, self.pose)

    @property
    def end_H4(self):
        return self.antibody_info.get_CDR_end(antibody.CDRNameEnum.H4, self.pose)

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
