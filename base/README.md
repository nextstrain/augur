## Original, non-modular augur

The older, original version of Augur is non-modular and doesn't use the `augur` command.
It is a set of Python 2 modules which can be used in scripts.
Currently, it's still lightly maintained to support older pathogen pipelines for Nextstrain that we have yet to move to the newest version of Augur.
However, we're actively migrating our pipelines to use the newest Augur and eventually old Augur will be removed.
For the time being, both versions live inside this git repository, which can be confusing.

These are the files associated with new, modular Augur:

    augur/
      __init__.py
      align.py
      export.py
      ...
    bin/
      augur
    setup.py

and these are the files associated with the old Augur:

    __init__.py
    base/
      auspice_export.py
      colorLogging.py
      config.py
      ...
      requirements.txt
      requirements-locked.txt
    builds/
      dengue/
      flu/
    scripts/
      beast-to-auspice-jsons-proof-of-principle.py
      json_tree_to_nexus.py
      plot_msa.py

Old Augur requires writing Python scripts which import functions and classes from the `base/` directory.
Examples of this are in some of our old pipelines in the `builds/` directory; see the `prepare.py` and `process.py` scripts.

Running old Augur requires Python 2.7, unlike the latest version of Augur, and dependencies are listed in `requirements.txt`.
You can install them with a package manager like conda or pip:

    pip install -r requirements.txt

You may choose to use `requirements-locked.txt` instead if you'd like a fixed set of known-good packages.
