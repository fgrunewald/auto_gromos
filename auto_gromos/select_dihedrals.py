from collections import defaultdict
from vermouth.processors.processor import Processor

class SelectDihedrals(Processor):
    """
    For each unique rotatable bond select
    one dihedral based on a priority scheme.
    Dihedrals with higher priority are more
    specific and chosen over dihedrals with lower
    priority. The priority is an attribute of the
    Interactions.meta dict. Note that improper
    dihedrals are excluded from this rule.
    """

    def run_molecule(self, molecule):
        """
        Select one dihedrals for each rotatable bond
        based on using the dihdral with higher priority.
        """
        if "dihedrals" in molecule.interactions:
            dih_dict = defaultdict(list)
            for dih in molecule.interactions["dihedrals"]:
                atoms = dih.atoms
                dih_dict[(frozenset([atoms[1], atoms[2]]))].append(dih)
            
            molecule.interactions["dihedrals"] = []
            
            for bond in dih_dict:
                prios = [dih.meta["priority"]
                         for idx, dih in enumerate(dih_dict[bond])]
                keep = prios.index(max(prios))
                molecule.interactions["dihedrals"].append(dih_dict[bond][keep])

        return molecule
