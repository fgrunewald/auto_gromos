import networkx as nx

def write_blocks(blocks, outname):
    """
    Write `blocks` in a file `outname`
    using vermouth format.

    Parameters:
    -----------
    blocks: list[:class:`vermouth.molecule.blocks`]
    outname: str
    """
    with open(outname, "w") as _file:
        for block in blocks:
           # write header
            _file.write("[ moleculetype ]\n")
            name = nx.get_node_attributes(block, "resname")[0]
            _file.write(name + " 3\n")
            _file.write("[ atoms ]\n")

            # write atom section
            for node in block.nodes:
                atom = block.nodes[node]
                _format_atoms = " ".join(['{}' for _ in atom.keys()]) + "\n"
                _line_values = [node+1] + [atom[key] for key in ["atype", "resid",
                                                               "resname", "atomname",
                                                               "charge_group", "charge", "mass"]]
                _file.write(_format_atoms.format(*_line_values))

            # write interactions and keep comments
            for inter_type in block.interactions:
                _file.write("[ " + inter_type + " ]\n")
                for inter in block.interactions[inter_type]:
                    _format_atoms = " ".join(['{}' for _ in inter.atoms])
                    _format_inter = " ".join(['{}' for _ in inter.parameters])

                    values = [atom+1 for atom in inter.atoms] + inter.parameters
                    if "comment" in inter.meta:
                        _format = _format_atoms + " " + _format_inter + \
                            " {{\"comment\":\"{}\"}}" + '\n'
                        values.append(inter.meta["comment"])
                    else:
                        _format = _format_atoms + " " + _format_inter + '\n'
                    _file.write(_format.format(*values))


def write_links(links, outname):
    """
    Write `links` in a file `outname`
    using vermouth format.

    Parameters:
    -----------
    links: list[:class:`vermouth.molecule.link`]
    outname: str
    """
    written_links = {}
    with open(outname, "w") as _file:

        for link in links:
            _file.write("[ link ]\n")
            _file.write("resname " + "\"" + link.name + "\"" + "\n")

            for inter_type in link.interactions:
                _file.write("[ " + inter_type + " ]\n")

                for inter in link.interactions[inter_type]:
                    _format_atoms = " ".join(['{}' for _ in inter.atoms])
                    _format_inter = " ".join(['{}' for _ in inter.parameters])

                    values = list(inter.atoms) + inter.parameters

                    if "comment" in inter.meta:
                        _format = _format_atoms + " " + _format_inter + \
                            " {{\"comment\":\"{}\"}}" + '\n'
                        values.append(inter.meta["comment"])
                    else:
                        _format = _format_atoms + " " + _format_inter + '\n'

                    _file.write(_format.format(*values))
                    written_links[inter.atoms] = inter.parameters
