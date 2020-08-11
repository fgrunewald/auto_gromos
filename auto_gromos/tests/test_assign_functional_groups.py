import pytest
from vermouth import molecule

@pytest.fixture
def test_molecule(functional_group):
    molecule = vermouth.molecule.Molecule()
    # The node keys should not be in a sorted order as it would mask any issue
    # due to the keys being accidentally sorted.
    molecule.add_node(2, atomname='CC')
    molecule.add_node(0, atomname='AA')
    molecule.add_node(1, atomname='BB')

    molecule.add_edge(0, 1)
    molecule.add_edge(0, 2)
    molecule.add_interaction(
        type_='bonds',
        atoms=(0, 1),
        parameters=['1', '2'],
        meta={'unmutable': 0, 'mutable': [4, 5, 6]},
    )
    molecule.add_interaction(
        type_='bonds',
        atoms=(0, 2),
        parameters=['a', 'b'],
    )
    return molecule
