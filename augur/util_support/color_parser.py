import functools
from collections import defaultdict
from typing import Dict, List, Optional, TextIO, Tuple

from augur.data import as_file
from augur.io.file import open_file
from augur.util_support.color_parser_line import ColorParserLine


class ColorParser:
    def __init__(self, *, mapping_filename: Optional[str], use_defaults: bool = True) -> None:
        self.mapping_filename = mapping_filename
        self.use_defaults = use_defaults

    @property
    @functools.lru_cache()
    def mapping(self) -> Dict[str, List[Tuple[str, str]]]:
        colors: Dict[str, List[Tuple[str, str]]] = {}

        if self.use_defaults:
            with as_file("colors.tsv") as file:
                with open_file(file) as defaults:
                    colors = {**colors, **self.parse_file(defaults)}

        if self.mapping_filename:
            with open_file(self.mapping_filename) as mapping:
                colors = {**colors, **self.parse_file(mapping)}

        return colors

    def parse_file(self, file: TextIO) -> Dict[str, List[Tuple[str, str]]]:
        file_mapping: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
        for pair in [ColorParserLine(line).pair() for line in file]:
            if pair is None:
                continue

            file_mapping[pair[0]].append(pair[1])

        return file_mapping
