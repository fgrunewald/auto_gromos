import networkx as nx
from vermouth.processors.do_links import DoLinks
from auto_gromos import AssignFunctionalGroups, SelectDihedrals, Gen14Pairs, RemoveDuplicates


def delete_atomname(molecule):
    """
    Delete the atomname attribute form a molecule
    and return it as dict.

    Parameters:
    -----------
    molecule: `:class:vermouth.molecule`

    Returns:
    --------
    molecule
        the molecule with no atomname attribute
    atomnames
        the atomname dict
    """
    atomnames = nx.get_node_attributes(molecule, "atomname")
    for node in molecule.nodes:
        del molecule.nodes[node]["atomname"]

    return molecule, atomnames


def gen_bonded_interactions(molecule):
    """
    Given a molecule, generate all bonded
    interactions according to the definition of 
    the force-field.
    """
    mol, atomnames = delete_atomname(molecule)

    nx.set_node_attributes(mol, {node: 1 for node in mol.nodes}, "resid")

    mol = AssignFunctionalGroups.run_molecule(mol)
    mol = DoLinks().run_molecule(mol)
    mol = SelectDihedrals().run_molecule(mol)
    mol = Gen14Pairs().run_molecule(mol)
    mol = RemoveDuplicates().run_molecule(mol)

    nx.set_node_attributes(mol, atomnames, "atomname")

    return molecule
