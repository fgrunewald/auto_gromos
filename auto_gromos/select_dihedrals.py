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
        dih_dict = defaultdict(list)
        if "dihedrals" in molecule.interactions:
            for dih in molecule.interactions["dihedrals"]:
                atoms = dih.atoms
                dih_dict[(frozenset([atoms[1], atoms[2]]))].append(("dihedrals", dih))

        if "impropers" in molecule.interactions:
            for dih in molecule.interactions["impropers"]:
                atoms = dih.atoms
                dih_dict[(frozenset([atoms[1], atoms[2]]))].append(("impropers", dih))

        molecule.interactions["dihedrals"] = []
        molecule.interactions["impropers"] = []

        for bond in dih_dict:
            special = False
            for inter_type, dih in dih_dict[bond]:
                if "keep" in dih.meta:
                    special=True
                    molecule.interactions[inter_type].append(dih)

            if not special:
                prios = [dih[1].meta["priority"]
                        for idx, dih in enumerate(dih_dict[bond])]
                keep = prios.index(max(prios))
                inter_type, dih_to_keep = dih_dict[bond][keep]
                molecule.interactions[inter_type].append(dih_to_keep)

        return molecule
