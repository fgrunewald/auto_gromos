import itertools
from collections import defaultdict
import vermouth
from vermouth.forcefield import ForceField
from vermouth.ffinput import read_ff
from vermouth.molecule import Interaction
import networkx as nx
from auto_gromos import DATA_PATH
from auto_gromos.assign_functional_groups import delete_atomname
from auto_gromos.gen_bonded_interactions import gen_bonded_interactions
from auto_gromos.ff_writers import write_links

def _treat_link_atoms(molecule, link, inter_type):
    # the uncommented statement does not work because node and
    # atom name are couple for blocks, which is debatably useful
    #atom_names = list(nx.get_node_attributes(block, 'atomname'))
    atom_names = nx.get_node_attributes(molecule, "atomname")
    resids = nx.get_node_attributes(molecule, "resid")
    new_interactions = defaultdict(list)
    for inter_type in link.interactions:
        for interaction in link.interactions[inter_type]:
            new_atoms = []
            for atom in interaction.atoms:
                order = resids[atom] - 1
                prefix = "".join(["+" for _ in range(0, order)])
                new_name = prefix + atom_names[atom]
                new_atoms.append(new_name)
                attrs = molecule.nodes[atom]
                link.add_node(new_name, **attrs)
                nx.set_node_attributes(link, {new_name: order}, "order")

            new_inter = Interaction(atoms=tuple(new_atoms),
                                    parameters=interaction.parameters,
                                    meta=interaction.meta)
            new_interactions[inter_type].append(new_inter)

    link.interactions.update(new_interactions)
    return new_atoms


def find_atoms(molecule, **attrs):
    """
    Yields all indices of atoms that match `attrs`
    Parameters
    ----------
    molecule: :class:`vermouth.molecule.Molecule`
    **attrs: collections.abc.Mapping
        The attributes and their desired values.
    Yields
    ------
    collections.abc.Hashable
        All atom indices that match the specified `attrs`
    """
    for node_idx in molecule:
        node = molecule.nodes[node_idx]
        if vermouth.molecule.attributes_match(node, attrs, ignore_keys=['resname', 'charge', 'charge_group', 'mass', 'atype']):
            yield node_idx


def extract_links(molecule, force_field):
    node_to_resid = nx.get_node_attributes(molecule, "resid")
    for inter_type in molecule.interactions:
        links = []
        prev_atoms = []
        for interaction in molecule.interactions[inter_type]:
            atoms = interaction.atoms
            resnames = set([molecule.nodes[node]["resname"] for node in atoms])
            if len(set([node_to_resid[node] for node in atoms])) > 1:
                if interaction.atoms != prev_atoms:
                    prev_atoms[:] = interaction.atoms
                    new_link = vermouth.molecule.Link()
                    new_link.interactions = defaultdict(list)
                    new_link.name = "|".join(list(resnames))
                    links.append(new_link)
                links[-1].interactions[inter_type].append(interaction)

        for link in links:
            _treat_link_atoms(molecule, link, inter_type)
            force_field.links.append(link)

    return force_field


def gen_links_gromos2016(force_field, names, graph_links, moltypes):
    """
    Generate all links for the `names` in `force_field`
    for the gromos2016 force_field. Note that graph
    links must contain definitions of links at graph
    level.
    """
    with open(DATA_PATH + "/2016H66/gromos2016H66_links.ff", "r") as _file:
        lines = _file.readlines()

    read_ff(lines, force_field)

    for link in force_field.links:
        delete_atomname(link)

    # generate combinations of dimers, trimers, etc.
    # this should be dynamic but even if there is only
    # 1 carbon atom in the backbone a tetramer captuers
    # the farthest interaction
    molecules = []
    for r in [1, 2, 3, 4]:
        for combo in itertools.combinations(names, r=r):
            if r == 1:
               combo = [combo[0], combo[0]]
            block_0 = force_field.blocks[combo[0]]
            prev_res = block_0.nodes[list(block_0.nodes())[0]]["resname"]
            prev_moltype = moltypes[prev_res]
            mol = block_0.to_molecule()
            resid = 1

            # apply all the edges between the blocks in the combination
            for name in combo[1:]:
                block = force_field.blocks[name]
                mol.merge_molecule(block)
                res = block.nodes[list(block.nodes())[0]]["resname"]
                moltype = moltypes[res]
                link = graph_links[(prev_moltype, moltype)]
                link[0].update({"resid": resid})
                link[1].update({"resid": resid+1})
                node1 = [i for i in find_atoms(mol, **link[0])][0]
                node2 = [i for i in find_atoms(mol, **link[1])][0]
                mol.add_edge(node1, node2)
                resid += 1
                prev_moltype = moltype

            # generate all bonded interactions for the fragment
            mol = gen_bonded_interactions(mol)
            molecules.append(mol)

    # extract links from the molecules
    polyply_links = ForceField("polyplylinks")
    for mol in molecules:
        extract_links(mol, polyply_links)

    # filter the links to not write duplicate links
    had_links = {}
    delete = []
    for idx, link in enumerate(polyply_links.links):
        unique_links = 0
        for inter_type in link.interactions:
            for inter in link.interactions[inter_type]:
                if (inter.atoms, link.name) in had_links:
                    if inter.parameters == had_links[(inter.atoms, link.name)][0]:
                        continue
                else:
                    unique_links += 1
                    had_links[(inter.atoms, link.name)] = (
                        inter.parameters, link.name)

        if unique_links == 0:
            delete.append(idx)

    delete.reverse()
    for idx in delete:
        del polyply_links.links[idx]

    # write links to file
    write_links(polyply_links.links, "gromos_links.ff")
