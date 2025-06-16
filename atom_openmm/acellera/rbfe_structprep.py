import sys


def rbfe_structprep(config_file=None):
    from atom_openmm.acellera.utils import set_directory, create_system
    from pathlib import Path
    import yaml
    import json
    import time
    from atom_openmm.rbfe_structprep import (
        massage_keywords,
        do_mintherm,
        do_lambda_annealing,
        do_equil,
    )
    import logging
    import os

    if config_file is None:
        config_file = sys.argv[1]

    print("")
    print("========================================")
    print("AToM RBFE Structure Preparation         ")
    print("========================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", config_file)
    print("")
    sys.stdout.flush()

    if config_file.endswith(".yaml"):
        with open(config_file, "r") as f:
            keywords = yaml.safe_load(f)
    elif config_file.endswith(".json"):
        with open(config_file, "r") as f:
            keywords = json.load(f)
    else:
        raise ValueError("Invalid configuration file format")

    logger = logging.getLogger("rbfe_structprep")

    with set_directory(Path(config_file).parent):
        create_system(
            os.path.join(keywords["BASENAME"] + ".prmtop"),
            os.path.join(keywords["BASENAME"] + ".inpcrd"),
            keywords,
        )
        restrain_solutes = True
        old_keywords = keywords.copy()
        massage_keywords(keywords, restrain_solutes)

        do_mintherm(keywords, logger)
        do_lambda_annealing(keywords, logger)

        # reestablish the restrained atoms
        if restrain_solutes:
            keywords["POS_RESTRAINED_ATOMS"] = old_keywords.get("POS_RESTRAINED_ATOMS")

        do_equil(keywords, logger)


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Specify ONE input file"

    rbfe_structprep(sys.argv[1])
