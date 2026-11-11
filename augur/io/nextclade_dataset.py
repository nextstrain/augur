"""
Source reference and annotation files from a Nextclade dataset.

A Nextclade dataset is a directory containing a ``pathogen.json`` whose
``files`` section names the other dataset files relative to the directory.
"""
import argparse
import json
import os
from dataclasses import dataclass
from typing import Optional

from augur.errors import AugurError


PATHOGEN_JSON = "pathogen.json"


@dataclass
class NextcladeDataset:
    path: str
    reference: str
    annotation: Optional[str]


def add_nextclade_dataset_argument(parser):
    """Add the --nextclade-dataset argument to *parser* (a parser, argument group or mutually exclusive group)."""
    return parser.add_argument('--nextclade-dataset', metavar="DIR", type=str,
        help="Nextclade dataset directory (containing a pathogen.json) from which to source the reference sequence and/or genome "
             "annotation instead of providing them individually. A genome annotation sourced from a dataset is read "
             "with --use-nextclade-gff-style.")


def read_nextclade_dataset(path: str) -> NextcladeDataset:
    """
    Read the pathogen.json of the Nextclade dataset at *path* (the dataset
    directory or the pathogen.json itself) and resolve the reference and
    genome annotation files.

    Raises
    ------
    AugurError
        If pathogen.json can't be read, doesn't declare a reference, or if a
        declared file doesn't exist
    """
    if os.path.isdir(path):
        dataset_dir = path
        pathogen_json = os.path.join(path, PATHOGEN_JSON)
    else:
        dataset_dir = os.path.dirname(path)
        pathogen_json = path

    if not os.path.isfile(pathogen_json):
        raise AugurError(f"Nextclade dataset {path!r} does not contain a {PATHOGEN_JSON}.")

    try:
        with open(pathogen_json, encoding="utf-8") as fh:
            pathogen = json.load(fh)
    except json.JSONDecodeError as error:
        raise AugurError(f"Could not parse {pathogen_json!r}: {error}")

    files = pathogen.get("files") if isinstance(pathogen, dict) else None
    if not isinstance(files, dict):
        raise AugurError(f"{pathogen_json!r} does not contain a 'files' section.")

    def resolve(key):
        if not files.get(key):
            return None
        file = os.path.join(dataset_dir, files[key])
        if not os.path.isfile(file):
            raise AugurError(f"File {files[key]!r} declared as files.{key} in {pathogen_json!r} does not exist.")
        return file

    reference = resolve("reference")
    if reference is None:
        raise AugurError(f"{pathogen_json!r} does not declare a reference sequence (files.reference).")

    return NextcladeDataset(path=path, reference=reference, annotation=resolve("genomeAnnotation"))


def apply_nextclade_dataset(args: argparse.Namespace, reference: Optional[str] = None, annotation: Optional[str] = None):
    """
    If ``args.nextclade_dataset`` is set, fill the argument with destination
    *reference* with the dataset's reference sequence and the argument with
    destination *annotation* with the dataset's genome annotation. Filling the
    annotation also sets ``args.use_nextclade_gff_style``.

    Raises
    ------
    AugurError
        If an argument to be filled was also provided explicitly, or if
        *annotation* is requested but the dataset has no genome annotation
    """
    if not args.nextclade_dataset:
        return
    dataset = read_nextclade_dataset(args.nextclade_dataset)

    for dest in (reference, annotation):
        if dest and getattr(args, dest):
            flag = "--" + dest.replace("_", "-")
            raise AugurError(f"--nextclade-dataset and {flag} can not be used together.")

    if reference:
        setattr(args, reference, dataset.reference)
    if annotation:
        if dataset.annotation is None:
            raise AugurError(f"Nextclade dataset {dataset.path!r} does not declare a genome annotation (files.genomeAnnotation).")
        setattr(args, annotation, dataset.annotation)
        args.use_nextclade_gff_style = True
