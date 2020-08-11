import networkx as nx
from vermouth.processors.processor import Processor
from vermouth.molecule import Choice, attributes_match
from vermouth.ffinput import read_ff
from vermouth.forcefield import ForceField
from auto_gromos import DATA_PATH

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

def _atoms_match(node1, node2):
    return attributes_match(node1, node2, ignore_keys=('order', 'replace', 'name'))


class AssignFunctionalGroups(Processor):
    """
    Given a molecule graph assign functinal groups
    which are needed for the generation of bonded
    interactions.
    """

    def _assign_functional_groups(self, func_groups):
        """
        Using the functinal group definitions of `func_groups`
        check if any of the subgraph is sub-isomorphic to the
        functional group graph and set the name of the functional
        group. Note that overlapping functional group definitions
        are allowed in order to facilitate for example assignment
        of tetrahedral centers.

        Parameters:
        -----------
        func_groups: `:class:vermouth.force_field.ForceField`
        """
        for group in func_groups.links:
            group, _ = delete_atomname(group)
            name = group.nodes['+0']["name"]
            GM = nx.isomorphism.GraphMatcher(
                self.molecule, group, node_match=_atoms_match)
            raw_matches = GM.subgraph_isomorphisms_iter()

            for raw_match in raw_matches:
                for node in raw_match.keys():
                    if "fgroup" in self.molecule.nodes[node]:
                        self.molecule.nodes[node].values.append(name)
                    else:
                        fgroup = Choice([name])
                        self.molecule.nodes[node]["fgroup"] = fgroup

    def _set_default_group(self):
        """
        If no functional group can be identified the fgroup
        attribute is set to None.
        """
        for node in self.molecule.nodes:
            if not "fgroup" in self.molecule.nodes[node]:
                self.molecule.nodes[node]["fgroup"] = Choice(["None"])

    def run_molecule(self, molecule):
        """
        Try to determine the functional group according to the
        definitions in the force-field file for each node/atom
        of the molecule. Whenever, no functional group is
        identified the fgroup attribute is set to None.
        """
        self.molecule = molecule
        func_groups = ForceField("functional-groups")
        # ToDo use pathlib for this
        with open(DATA_PATH + "gromos2016H66_fgroups.ff", "r") as _file:
            lines = _file.readlines()

        read_ff(lines, func_groups)

        self._assign_functional_groups(func_groups)
        self._set_default_group()

        return molecule
