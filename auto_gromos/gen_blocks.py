from vermouth.ffinput import read_ff
from auto_gromos import DATA_PATH
from auto_gromos.ff_writers import write_blocks
from auto_gromos.assign_functional_groups import delete_atomname
from auto_gromos.gen_bonded_interactions import gen_bonded_interactions

def gen_block_gromos2016(force_field, names):
    """
    Given definitions of blocks in `force_field`
    generate all missing bonded interactions
    according to the gromos 2016H66 force-field for
    the molecules provided in `names`.
    """

    # parked this here to not have so many force_field imports.
    # can be moved though
    with open(DATA_PATH + "/2016H66/gromos2016H66_links.ff", "r") as _file:
        lines = _file.readlines()

    read_ff(lines, force_field)

    for link in force_field.links:
        delete_atomname(link)

    molecules = []
    for name in names:
        block = force_field.blocks[name]
        mol = block.to_molecule()
        mol = gen_bonded_interactions(mol)
        molecules.append(mol)

    write_blocks(molecules, "gromos2016_blocks.ff")
