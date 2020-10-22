from collections import OrderedDict
import pytest
from vermouth.gmx import read_rtp
from vermouth.ffinput import read_ff
from vermouth.forcefield import ForceField
from vermouth.molecule import Molecule
from auto_gromos.gen_bonded_interactions import gen_bonded_interactions
from auto_gromos.assign_functional_groups import delete_atomname
from auto_gromos import DATA_PATH

def ref_molecules():

    with open("test_set.rtp") as _file:
        lines = _file.readlines()

    force_field = ForceField("ref-mols")
    read_rtp(lines, force_field)

    return force_field


def load_links():

    with open(DATA_PATH + "/2016H66/gromos2016H66_links.ff", "r") as _file:
        lines = _file.readlines()

    force_field = ForceField("links")
    read_ff(lines, force_field)

    for link in force_field.links:
        delete_atomname(link)

    return force_field

def _interaction_equal(interaction1, interaction2, inter_type):
    """
    Returns True if interaction1 == interaction2, ignoring rounding errors in
    interaction parameters.
    """
    p1 = list(map(str, interaction1.parameters))
    p2 = list(map(str, interaction2.parameters))
    a1 = list(interaction1.atoms)
    a2 = list(interaction2.atoms)
    if inter_type in ["bonds", "exclusions", "pairs"]:
        a1.sort()
        a2.sort()
        return a1 == a2 and p1 == p2

    elif inter_type in ["impropers", "dihedrals"]:
         if frozenset([a1[0], a1[1], a1[2]]) == frozenset([a2[0], a2[1], a2[2]]) and a1[3] == a2[3] and p1 == p2:
            return True
         elif frozenset([a1[0], a1[1], a1[2]]) == frozenset([a2[3], a2[2], a2[1]]) and a1[3] == a2[0] and p1 == p2:
            return True
         elif frozenset([a1[1], a1[2], a1[3]]) == frozenset([a2[1], a2[2], a2[3]]) and a1[0] == a2[0] and p1 == p2:
            return True

    elif inter_type in ["angles"]:
        return a1[1] == a2[1] and frozenset([a1[0], a1[2]]) == frozenset([a2[0], a2[2]]) and p1 == p2

    return False
@pytest.mark.parametrize("name, block", ref_molecules().blocks.items())
def test_gen_bonded_interactions(name, block):
    molecule = block.to_molecule()
    print(block.interactions.keys())
    new_mol = Molecule(force_field=load_links(), nrexcl=3)
    new_mol.name = name
    new_mol.add_nodes_from(molecule.nodes(data=True))
    new_mol.add_edges_from(molecule.edges)
    new_mol = gen_bonded_interactions(new_mol)
    print(new_mol.interactions.keys())

    for inter_type in ["bonds", "angles", "constraints", "exclusions", "dihedrals", "impropers"]:
        print(inter_type)
        interactions = molecule.interactions.get(inter_type, [])
        new_interactions = new_mol.interactions.get(inter_type, [])
        #print(interactions)
        #print(new_interactions)
        for inter in interactions:
            print(inter)
        for inter in new_interactions:
            print(inter)
        #print("ref", interactions, len(interactions))
        #print("new", new_interactions, len(new_interactions))
        assert len(interactions) == len(new_interactions)
        for inter in interactions:
            atoms = inter.atoms
            for inter_new in new_interactions:
                if set(atoms) == set(inter_new.atoms):
                   assert _interaction_equal(inter, inter_new, inter_type)

    #assert False
