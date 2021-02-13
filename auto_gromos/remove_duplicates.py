from collections import defaultdict
from vermouth.processors.processor import Processor

class RemoveDuplicates(Processor):
    """
    For all interactions in a molecule, check if an
    interaction is identified more than once. For
    example a bond 1-2 and a bond 2-1. In this case
    remove one of the redundant interactions.
    """

    def run_molecule(self, molecule):
        """
        Check for each interaction if it is already
        defined. If it is already defined keep the first
        interaction in the list.
        """
        for inter_type in molecule.interactions:
            inter_dict = defaultdict(list)
            for inter in molecule.interactions[inter_type]:
                if "version" in inter.meta:
                    version = inter.meta["version"]
                    inter_dict[(frozenset(inter.atoms), version)] = inter
                else:
                    inter_dict[(frozenset(inter.atoms), 1)] = inter

            molecule.interactions[inter_type] = []
            for inter in inter_dict.values():
                molecule.interactions[inter_type].append(inter)

        return molecule
