import shutil
import os


curr_dir = os.path.dirname(os.path.abspath(__file__))


def _test_abfe_structprep(tmp_path):
    from atom_openmm.abfe_structprep import abfe_structprep

    run_dir = os.path.join(tmp_path, "1STP")
    shutil.copytree(os.path.join(curr_dir, "1STP"), run_dir)

    abfe_structprep(os.path.join(run_dir, "1STP_asyncre.yaml"))


def _test_abfe_production(tmp_path):
    from atom_openmm.abfe_production import abfe_production

    shutil.copytree(
        os.path.join(curr_dir, "1STP_equil"),
        os.path.join(tmp_path, "1STP_equil"),
    )
    run_dir = os.path.join(tmp_path, "1STP_equil")
    abfe_production(os.path.join(run_dir, "1STP_asyncre.yaml"))


def _test_rbfe_structprep(tmp_path):
    from atom_openmm.rbfe_structprep import rbfe_structprep
    from openmm import unit
    from openmm import app
    from openmm.app.amberprmtopfile import AmberPrmtopFile
    from openmm import XmlSerializer
    import numpy as np

    run_dir = os.path.join(tmp_path, "QB_A08_A07")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07"), run_dir)

    box_vectors = np.array([[80.1356, 0, 0], [0, 86.0443, 0], [0, 0, 81.4926]])
    a = unit.Quantity(box_vectors[0] * unit.angstrom)
    b = unit.Quantity(box_vectors[1] * unit.angstrom)
    c = unit.Quantity(box_vectors[2] * unit.angstrom)
    box_vectors = (a, b, c)
    omm_structure = AmberPrmtopFile(
        os.path.join(run_dir, "QB_A08_A07.prmtop"), periodicBoxVectors=box_vectors
    )
    system = omm_structure.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=9 * unit.angstrom,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
        rigidWater=True,
        switchDistance=7 * unit.angstrom,
    )
    with open(os.path.join(run_dir, "QB_A08_A07_sys.xml"), "w") as f:
        f.write(XmlSerializer.serialize(system))

    rbfe_structprep(os.path.join(run_dir, "QB_A08_A07_asyncre.yaml"))


def _test_rbfe_production(tmp_path):
    from atom_openmm.rbfe_production import rbfe_production

    shutil.copytree(
        os.path.join(curr_dir, "QB_A08_A07_equil"),
        os.path.join(tmp_path, "QB_A08_A07_equil"),
    )
    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil")
    rbfe_production(os.path.join(run_dir, "QB_A08_A07_asyncre.yaml"))


def _test_make_atm_system_from_amber(tmp_path):
    from atom_openmm.make_atm_system_from_amber import make_system

    refxml = os.path.join(curr_dir, "QB_A08_A07_equil", "QB_A08_A07_sys.xml")
    outxml = os.path.join(tmp_path, "QB_A08_A07_sys.xml")
    prmtopfile = os.path.join(curr_dir, "QB_A08_A07", "QB_A08_A07.prmtop")
    crdfile = os.path.join(curr_dir, "QB_A08_A07", "QB_A08_A07.inpcrd")
    pdboutfile = os.path.join(tmp_path, "QB_A08_A07_sys.pdb")
    make_system(
        prmtopfile=prmtopfile,
        crdfile=crdfile,
        xmloutfile=outxml,
        pdboutfile=pdboutfile,
        hmass=1.5,
        switchDistance=0.7,
        nonbondedCutoff=0.9,
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    with open(refxml, "r") as f:
        ref_lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"

    os.remove(outxml)
    os.system(
        f"make_atm_system_from_amber --AmberPrmtopinFile {prmtopfile} --AmberInpcrdinFile {crdfile} --systemXMLoutFile {outxml} --systemPDBoutFile {pdboutfile} --hmass 1.5 --switchDistance 0.7 --nonbondedCutoff 0.9"
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"


def _test_make_atm_system_from_pdb(tmp_path):
    from atom_openmm.make_atm_system_from_pdb import make_system

    refxml = os.path.join(curr_dir, "3ptb", "3ptb_sys.xml")
    outxml = os.path.join(tmp_path, "3ptb_sys.xml")
    pdb = os.path.join(curr_dir, "3ptb", "3ptb.pdb")
    sdffile = os.path.join(curr_dir, "3ptb", "BEN_ideal.sdf")
    pdboutfile = os.path.join(tmp_path, "3ptb_sys.pdb")
    make_system(
        systempdbfile=pdb,
        ligandsdffile=sdffile,
        lig1resid=1,
        xmloutfile=outxml,
        pdboutfile=pdboutfile,
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    with open(refxml, "r") as f:
        ref_lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"

    os.remove(outxml)
    os.system(
        f"make_atm_system_from_pdb --systemPDBinFile {pdb} --ligandsSDFFile {sdffile} --LIG1resid 1 --systemXMLoutFile {outxml} --systemPDBoutFile {pdboutfile}"
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"


def _test_make_atm_system_from_rcpt_lig(tmp_path):
    from atom_openmm.make_atm_system_from_rcpt_lig import make_system

    refxml = os.path.join(curr_dir, "3ptb", "3ptb_sys_2.xml")
    outxml = os.path.join(tmp_path, "3ptb_sys_2.xml")
    pdb = os.path.join(curr_dir, "3ptb", "3ptb.pdb")
    sdffile = os.path.join(curr_dir, "3ptb", "BEN_ideal.sdf")
    pdboutfile = os.path.join(tmp_path, "3ptb_sys_2.pdb")
    make_system(
        receptorfile=pdb,
        lig1sdffile=sdffile,
        displacement=[22, 22, 22],
        xmloutfile=outxml,
        pdboutfile=pdboutfile,
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    with open(refxml, "r") as f:
        ref_lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"

    os.remove(outxml)
    os.system(
        f"make_atm_system_from_rcpt_lig --receptorinFile {pdb} --LIG1SDFinFile {sdffile} --displacement '22.0 22.0 22.0' --systemXMLoutFile {outxml} --systemPDBoutFile {pdboutfile}"
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"


def _test_input_parser():
    from atom_openmm.utils.config import parse_config

    config_yaml = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.yaml"))
    config_json = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.json"))
    config_cntl = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.cntl"))

    assert config_yaml == config_json
    assert config_yaml == config_cntl

    for key in config_yaml.keys():
        assert config_yaml[key] == config_json[key]
        assert config_yaml[key] == config_cntl[key]
        assert type(config_yaml[key]) is type(
            config_json[key]
        ), f"{key} {type(config_yaml[key])} {type(config_json[key])}"
        assert type(config_yaml[key]) is type(
            config_cntl[key]
        ), f"{key} {type(config_yaml[key])} {type(config_cntl[key])}"
