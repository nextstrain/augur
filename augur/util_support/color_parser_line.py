import functools
import re


class ColorParserLine:
    def __init__(self, line):
        self.line = line

    def pair(self):
        if self.is_comment_or_blank():
            return None

        if len(self.fields) != 3:
            print("WARNING: Color map file contains invalid line:", self.line)
            return None

        if not self.hex_code.startswith("#") or len(self.hex_code) != 7:
            print(
                "WARNING: Color map file contained this invalid hex code: ",
                self.hex_code,
            )
            return None

        return self.trait, (self.trait_value, self.hex_code)

    def is_comment_or_blank(self):
        return bool(self.line.startswith("#") or re.match(r"^\s*$", self.line))

    @property
    @functools.lru_cache()
    def fields(self):
        return self.line.strip().split("\t")

    @property
    def trait(self):
        return self.fields[0].lower()

    @property
    def trait_value(self):
        return self.fields[1].lower()

    @property
    def hex_code(self):
        return self.fields[2]
