import networkx as nx
from vermouth.processors.do_links import DoLinks
from auto_gromos import AssignFunctionalGroups, SelectDihedrals, Gen14Pairs, RemoveDuplicates
from auto_gromos.assign_functional_groups import delete_atomname

def gen_bonded_interactions(molecule):
    """
    Given a molecule, generate all bonded
    interactions according to the definition of 
    the force-field.
    """
    mol, atomnames = delete_atomname(molecule)

    resids = nx.get_node_attributes(mol, "resid")
    nx.set_node_attributes(mol, {node: 1 for node in mol.nodes}, "resid")

    mol = AssignFunctionalGroups().run_molecule(mol)
    mol = DoLinks().run_molecule(mol)
    mol = SelectDihedrals().run_molecule(mol)
    mol = Gen14Pairs().run_molecule(mol)
    mol = RemoveDuplicates().run_molecule(mol)

    nx.set_node_attributes(mol, atomnames, "atomname")
    nx.set_node_attributes(mol, resids, "resid")

    return molecule
