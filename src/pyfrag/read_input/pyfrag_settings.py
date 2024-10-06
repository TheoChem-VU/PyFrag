import pathlib as pl
from typing import List, Protocol, Union

from pydantic import BaseModel, Field


class PyFragInputKey(Protocol):
    def label(self) -> str: ...


class BondLength(BaseModel):
    atom1: int
    atom2: int
    bond_length: float

    def label(self) -> str:
        return f"{self.atom1}-{self.atom2}"


class Angle(BaseModel):
    atom1: int
    atom2: int
    angle: float

    def label(self) -> str:
        return f"{self.atom1}-{self.atom2}"


class Dihedral(BaseModel):
    atom1: int
    atom2: int
    atom3: int
    dihedral_angle: float

    def label(self) -> str:
        return f"{self.atom1}-{self.atom2}-{self.atom3}"


class Overlap(BaseModel):
    frag1: str
    MO1: str
    frag2: str
    MO2: str

    def label(self) -> str:
        return f"{self.frag1}-{self.MO1}-{self.frag2}-{self.MO2}"


class OverlapWithIrrep(BaseModel):
    irrep1: str
    frag1: str
    index1: str
    irrep2: str
    frag2: str
    index2: str

    def label(self) -> str:
        return f"{self.irrep1}-{self.frag1}-{self.index1}-{self.irrep2}-{self.frag2}-{self.index2}"


class Population(BaseModel):
    frag1: str
    MO1: str

    def label(self) -> str:
        return f"{self.frag1}-{self.MO1}"


class PopulationWithIrrep(BaseModel):
    irrep1: str
    frag1: str
    index1: str

    def label(self) -> str:
        return f"{self.irrep1}-{self.frag1}-{self.index1}"


class OrbitalEnergy(BaseModel):
    frag1: str
    MO1: str

    def label(self) -> str:
        return f"{self.frag1}-{self.MO1}"


class OrbitalEnergyWithIrrep(BaseModel):
    irrep1: str
    frag1: str
    index1: str

    def label(self) -> str:
        return f"{self.irrep1}-{self.frag1}-{self.index1}"


class PyFragSection(BaseModel):
    bondlengths: List[BondLength] = Field(default_factory=list)
    angles: List[Angle] = Field(default_factory=list)
    dihedrals: List[Dihedral] = Field(default_factory=list)
    overlaps: List[Union[Overlap, OverlapWithIrrep]] = Field(default_factory=list)
    populations: List[Union[Population, PopulationWithIrrep]] = Field(default_factory=list)
    orbitalenergies: List[Union[OrbitalEnergy, OrbitalEnergyWithIrrep]] = Field(default_factory=list)
    vdd_indices: List[int] = Field(default_factory=list)
    strain_values: List[float] = Field(default_factory=list)
    fragment_indices: List[List[int]] = Field(default_factory=list)
    coordfile: List[pl.Path] = Field(default_factory=list)
