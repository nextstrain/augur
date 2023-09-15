from collections import defaultdict
from importlib.resources import open_text
import functools

from augur.util_support.color_parser_line import ColorParserLine


class ColorParser:
    def __init__(self, *, mapping_filename, use_defaults=True):
        self.mapping_filename = mapping_filename
        self.use_defaults = use_defaults

    @property
    @functools.lru_cache()
    def mapping(self):
        colors = {}

        if self.use_defaults:
            with open_text("augur.data", "colors.tsv") as defaults:
                    colors = {**colors, **self.parse_file(defaults)}

        if self.mapping_filename:
            with open(self.mapping_filename, encoding="utf-8") as mapping:
                colors = {**colors, **self.parse_file(mapping)}

        return colors

    def parse_file(self, file):
        file_mapping = defaultdict(list)
        for pair in [ColorParserLine(line).pair() for line in file]:
            if pair is None:
                continue

            file_mapping[pair[0]].append(pair[1])

        return file_mapping
