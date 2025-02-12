"""
Export JSON files suitable for visualization with Auspice.
The JSON schema is available at <https://nextstrain.org/schemas/dataset/v2>
"""
# FIXME: re-exporting these from export.py to avoid rewiring top-level augur
# command and show a useful diff. I think it'd make more sense to just remove
# this file and rename export_v2.py to export.py (though this may not show up as
# a rename in the squashed PR diff, but renaming in a commit should preserve
# history)
from .export_v2 import register_parser, run
