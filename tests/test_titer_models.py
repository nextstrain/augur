from augur.reconstruct_sequences import load_alignments
from augur.titer_model import TiterCollection, TreeModel, SubstitutionModel
from augur.utils import read_tree


def test_titer_collection():
    # Confirm that titers load from a file path.
    titers = TiterCollection("tests/data/titer_model/h3n2_titers_subset.tsv")

    # Confirm that all distinct test and reference strains have been counted.
    assert len(titers.strains) == 62


def test_titer_tree_model_validate_measurements():
    """Validate titer tree model with a subset of measurements.
    """
    # Load tree.
    tree = read_tree("tests/data/titer_model/h3n2_ha_tree.nwk")

    # Prepare model.
    model = TreeModel(tree, "tests/data/titer_model/h3n2_titers_subset.tsv")
    model.prepare(
        training_fraction=0.8,
        subset_strains=False,
    )

    # Train the model.
    model.train()

    # Validate model on a subset of measurements.
    performance = model.validate()


def test_titer_tree_model_validate_strains():
    """Validate titer tree model with a subset of virus strains.
    """
    # Load tree.
    tree = read_tree("tests/data/titer_model/h3n2_ha_tree.nwk")

    # Prepare model.
    model = TreeModel(tree, "tests/data/titer_model/h3n2_titers_subset.tsv")
    model.prepare(
        training_fraction=0.8,
        subset_strains=True,
    )

    # Train the model.
    model.train()

    # Validate model on a subset of measurements.
    performance = model.validate()


def test_titer_substitution_model_validate_measurements():
    """Validate titer substitution model with a subset of measurements.
    """
    # Load alignments.
    alignments = load_alignments(["tests/data/titer_model/h3n2_ha_aligned_genbank_HA1.fasta"], ["HA1"])

    # Prepare model.
    model = SubstitutionModel(alignments, "tests/data/titer_model/h3n2_titers_subset.tsv")
    model.prepare(
        training_fraction=0.8,
        subset_strains=False,
    )

    # Train the model.
    model.train()

    # Validate model on a subset of measurements.
    performance = model.validate()
