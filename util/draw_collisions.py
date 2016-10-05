"""
Created on 27/04/2013

@author: thom
"""

import math

from rdkit.Chem import AllChem as Chem

from atoms.kinetic_molecule import KineticMolecule
from reactions.chemistry_factory import ChemistryFactory
from reactions.emergent_reactions import EmergentReactions
from reactors.spatial_reactor import SpatialReactor


def draw(reactant_mols, product_mols, label=""):

    in_v = [mol.get_velocity() for mol in reactant_mols]
    out_v = [mol.get_velocity() for mol in product_mols]

    from itertools import chain
    c = list(chain.from_iterable(in_v))
    c.extend(chain.from_iterable(out_v))
    max_scale = max(c)
    min_scale = min(c)
    scale = max(max_scale, math.copysign(min_scale, max_scale))

    in_v_radial = [vessel._xyz_to_radial(*v) for v in in_v]
    out_v_radial = [vessel._xyz_to_radial(*v) for v in out_v]

    caption = " + ".join(["{}".format(Chem.MolToSmiles(Chem.RemoveHs(mol))) for mol in reactant_mols])
    caption += r"$\rightarrow$ " + " + ".join(["{}".format(Chem.MolToSmiles(Chem.RemoveHs(mol))) for mol in product_mols])

    tikz = r'\subcaptionbox{{{}}}'.format(caption) + r'[0.5\textwidth]{%'
    tikz += r"""
\begin{tikzpicture}[scale=3.75,tdplot_main_coords]
\coordinate (O) at (0,0,0);
\draw[thick,->] (0,0,0) -- (1,0,0) node[anchor=north east]{$x$};
\draw[thick,->] (0,0,0) -- (0,1,0) node[anchor=north west]{$y$};
\draw[thick,->] (0,0,0) -- (0,0,1) node[anchor=south]{$z$};"""
    tikz += "\n" + r"\node at (1,1,0) {{{}}};".format(label) + "\n"

    initial_CM_radial = vessel._xyz_to_radial(*vessel._get_CM_velocity(reactant_mols))  # azimuth, pi/4-polar, r
    CM_radial = initial_CM_radial

    tikz += r'\tdplotsetcoord{{CM0}}{{{2}}}{{{1}}}{{{0}}};\tdplotsetcoord{{CM1}}{{{3}}}{{{1}}}{{{0}}};' .format(
        str(CM_radial[0] * 180.0), str(CM_radial[1] * 180.0), str(CM_radial[2] / scale), str(-CM_radial[2] / scale))
    tikz += '\n' + r'\draw[color=black,densely dotted,-stealth] (CM0)-- (CM1);'

    I_point_template = r'\tdplotsetcoord{{I{0}}}{{{3}}}{{{2}}}{{{1}}};'  # r, polar, azimuth
    I_draw_template = r'\draw[color=red,-stealth] (I{}) node[font=\tiny,color=black]{{{}}}-- (O);'
    I_projection_template = r'\draw[dotted, color=black] (O) -- (I{0}xy); \draw[dotted, color=black] (I{0}) -- (I{0}xy);'

    tikz += "\n".join([I_point_template.format(str(n), str(v[0] * 180.0), str(v[1] * 180.0), str(-v[2] / scale))
                       for n, v in zip(range(len(in_v_radial)), in_v_radial)]) + "\n"  # NEGATIVE r as coming FROM In to Origin
    tikz += "\n".join([I_draw_template.format(str(n), Chem.MolToSmiles(mol)) for n, mol in zip(range(len(in_v_radial)), reactant_mols)]) + "\n"
    tikz += "\n".join([I_projection_template.format(str(n)) for n in range(len(in_v_radial))]) + "\n"

    O_point_template = r'\tdplotsetcoord{{O{0}}}{{{3}}}{{{2}}}{{{1}}};'
    O_draw_template = r'\draw[color=blue,-stealth] (O) -- (O{}) node[font=\tiny,color=black]{{{}}};'
    O_projection_template = r'\draw[dotted, color=black] (O) -- (O{0}xy); \draw[dotted, color=black] (O{0}) -- (O{0}xy);'

    tikz += "\n".join([O_point_template.format(str(n), str(v[0] * 180.0), str(v[1] * 180.0), str(v[2] / scale))
                       for n, v in zip(range(len(out_v_radial)), out_v_radial)]) + "\n"  # POSITIVE r as TO On from Origin
    tikz += "\n".join([O_draw_template.format(str(n), Chem.MolToSmiles(mol)) for n, mol in zip(range(len(out_v_radial)), product_mols)]) + "\n"
    tikz += "\n".join([O_projection_template.format(str(n)) for n in range(len(out_v_radial))]) + "\n"

    tikz += r'\end{tikzpicture}}%' + '\n'

    return tikz


if __name__ == "__main__":

    vessel = SpatialReactor(-7, dimension=3)
    vessel.initialize(chemistry=ChemistryFactory.new())

    with open("Collisions.tex", "w") as f:
        f.write(r"""\tdplotsetmaincoords{60}{110}
\begin{figure}
\centering
""")

        reactant_mols = [KineticMolecule('O=C=O', kinetic_energy=300), KineticMolecule('C', kinetic_energy=200)]
        options = EmergentReactions(ChemistryFactory.new()).get_reaction_options(reactant_mols)
        diagram = 1
        for rxn in options:
            product_mols = rxn.fire()

            out_v, IE = vessel._inelastic_collision(reactant_mols, product_mols, 20)
            for mol, v in zip(product_mols, out_v):
                mol.set_velocity(*v)

            f.write(draw(reactant_mols, product_mols))
            if diagram % 2 == 0 and diagram < len(options):
                f.write(r"""\end{figure}
\begin{figure}
\centering
""")
            diagram += 1

        f.write(r'\end{figure}')
