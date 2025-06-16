from contextlib import contextmanager
from pathlib import Path
import os


@contextmanager
def set_directory(path: Path):
    """Sets the cwd within the context

    Args:
        path (Path): The path to the cwd

    Yields:
        None
    """

    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def create_system(prmtopfile, crdfile, keywords):
    from openmm.unit import amu, nanometer
    from openmm import app, XmlSerializer
    from moleculekit.molecule import Molecule
    import os
    import logging

    logger = logging.getLogger("atom_openmm.ommsystem")

    mol = Molecule(prmtopfile)
    mol.read(crdfile)
    mol.write(os.path.basename(prmtopfile).replace(".prmtop", ".pdb"))

    prmtop = app.AmberPrmtopFile(prmtopfile)
    if keywords.get("HMASS") is not None:
        hmass = float(keywords.get("HMASS")) * amu
    else:
        hmass = 1.0 * amu
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * nanometer,
        constraints=app.HBonds,
        hydrogenMass=hmass,
    )
    topology = prmtop.topology

    if nnp_model := keywords.get("NNP_MODEL"):
        from openmmml import MLPotential

        logger.info("Initialize NNP/MM")

        group_indices = []
        for i in [1, 2]:
            atom_indices = list(map(int, keywords[f"LIGAND{i}_ATOMS"]))
            logger.info(f"Ligand{i}: NNP atom indices: {atom_indices}")
            logger.info(f"Ligand{i}: number of NNP atoms: {len(atom_indices)}")
            assert len(atom_indices) > 0
            group_indices.append(atom_indices)

        logger.info(f"NNP model: {nnp_model}")
        nnp = None
        if nnp_model.startswith("TorchMD-NET"):
            from atom_openmm.acellera import atom_nnp_wrapper

            nnp_file = keywords["NNP_FILE"]
            logger.info(f"NNP file: {nnp_file}")
            max_num_neighbors = keywords["NNP_MAX_NUM_NEIGHBORS"]
            logger.info(f"NNP max num neighbors: {max_num_neighbors}")
            nnp = MLPotential(
                nnp_model,
                model_file=nnp_file,
                group_indices=group_indices,
                max_num_neighbors=max_num_neighbors,
                use_cuda_graphs=True,
            )
        else:
            nnp = MLPotential(nnp_model)

        all_atom_indices = group_indices[0] + group_indices[1]
        system = nnp.createMixedSystem(
            topology, system, all_atom_indices, removeConstraints=False
        )

    with open(os.path.basename(prmtopfile).replace(".prmtop", "_sys.xml"), "w") as f:
        f.write(XmlSerializer.serialize(system))
