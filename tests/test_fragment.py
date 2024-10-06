import pytest

# Import the necessary classes and functions
from pyfrag.jobs.fragment import Complex, Fragment, create_fragments_and_complex
from pyfrag.read_input.inputblocks import InputBlocks
from pyfrag.read_input.pyfrag_settings import PyFragSection
from scm.plams import Molecule, Settings


@pytest.fixture
def complex_instance():
    return Complex(name="test_complex", fragments=[], trajectory=[Molecule()], plams_settings=Settings())


@pytest.fixture
def fragment_instance():
    return Fragment(name="test_fragment", index=0, strain_energy=0.0, trajectory=[Molecule()], plams_settings=Settings())


def test_add_fragment(complex_instance, fragment_instance):
    complex_instance.add_fragment(fragment_instance)
    assert len(complex_instance.fragments) == 1
    assert complex_instance.fragments[0] == fragment_instance


def test_get_fragment(complex_instance, fragment_instance):
    complex_instance.add_fragment(fragment_instance)
    retrieved_fragment = complex_instance.get_fragment(0)
    assert retrieved_fragment == fragment_instance


def test_create_fragments_and_complex():
    file_blocks = InputBlocks(FRAGMENT_EXTRA={1: "", 2: ""})
    complex_name = "test_complex"
    frag_names = ["frag1", "frag2"]
    pyfrag_keys = PyFragSection(strain_values=[0.0, 1.0])

    complex, fragments = create_fragments_and_complex(
        file_blocks=file_blocks, mol_trajectories=[[Molecule()], [Molecule()], [Molecule()]], complex_name=complex_name, frag_names=frag_names, pyfrag_keys=pyfrag_keys
    )

    assert complex.name == complex_name
    assert len(fragments) == 2
    assert fragments[0].name == frag_names[0]
    assert fragments[1].name == frag_names[1]
    assert fragments[0].strain_energy == 0.0
    assert fragments[1].strain_energy == 1.0


if __name__ == "__main__":
    pytest.main()
