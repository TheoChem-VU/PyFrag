from typing import List, Tuple, Union

from pydantic import BaseModel, ConfigDict, Field
from scm.plams import Molecule, Settings

from pyfrag.process_input.process_ams import convert_input_block_to_settings, process_ams_input
from pyfrag.read_input.inputblocks import InputBlocks
from pyfrag.read_input.pyfrag_settings import PyFragSection


class Complex(BaseModel):
    name: str
    trajectory: List[Molecule]
    plams_settings: Settings
    fragments: List["Fragment"] = Field(default_factory=list)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def add_fragment(self, fragment: "Fragment") -> None:
        self.fragments.append(fragment)

    def get_fragment(self, index: int) -> "Fragment":
        return self.fragments[index]


class Fragment(BaseModel):
    name: str
    index: int
    strain_energy: float
    trajectory: List[Molecule]
    plams_settings: Settings

    model_config = ConfigDict(arbitrary_types_allowed=True)


def create_fragments_and_complex(
    file_blocks: InputBlocks, mol_trajectories: List[List[Molecule]], complex_name: Union[str, None], frag_names: Union[List[str], None], pyfrag_keys: PyFragSection
) -> Tuple[Complex, List[Fragment]]:
    """Function that creates the fragments and complex from the input file."""

    common_settings, calc_type = process_ams_input(file_blocks)

    # First, create the complex
    complex_name = "full" if complex_name is None else complex_name
    complex_traj = mol_trajectories[0]
    complex_settings = convert_input_block_to_settings(file_blocks.COMPLEX_EXTRA)
    complex = Complex(name=complex_name, trajectory=complex_traj, plams_settings=common_settings + complex_settings)

    # Then, create the fragments
    fragments = []
    frag_names = [f"f{i+1}" for i in range(len(mol_trajectories))] if frag_names is None else frag_names

    for i in range(len(mol_trajectories) - 1):
        frag_name = frag_names[i]
        strain_energy = pyfrag_keys.strain_values[i] if i > 0 else 0.0
        plams_settings = convert_input_block_to_settings(file_blocks.FRAGMENT_EXTRA.get(i + 1, ""))

        fragment = Fragment(name=frag_name, index=i, strain_energy=strain_energy, trajectory=mol_trajectories[i + 1], plams_settings=common_settings + plams_settings)
        complex.add_fragment(fragment)
        fragments.append(fragment)

    return complex, fragments
