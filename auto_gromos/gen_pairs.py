import networkx as nx
from vermouth.processors.processor import Processor
from vermouth.molecule import Interaction

def find_neighbours(molecule, source, max_length, min_length):
    paths = nx.single_source_shortest_path(G=molecule, source=source, cutoff=max_length)
    neighbours = [node for node, path in paths.items() if min_length <= len(path)]
    return neighbours

class Gen14Pairs(Processor):
    """
    Generate all 1-4 pairs for a molecule
    and add them to the interaction dict.
    """

    def run_molecule(self, molecule):
        """
        Given the graph of molecule find all
        1-4 pairs and add them to the interaction
        list.
        """
        for node in molecule.nodes:
            pairs = find_neighbours(molecule, node, 4, 4)
            for ndx in pairs:
                molecule.interactions["pairs"].append(
                    Interaction(atoms=(ndx, node), parameters=["2"], meta={}))
