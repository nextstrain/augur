from augur.filtering.matchers.date import Date
from augur.filtering.matchers.length import Length
from augur.filtering.matchers.metadata import Metadata
from augur.filtering.matchers.name import Name
from augur.filtering.matchers.nonnucleotide import Nonnucleotide


MATCHER_CLASSES = {
    "date": Date,
    "length": Length,
    "metadata": Metadata,
    "name": Name,
    "non-nucleotide": Nonnucleotide,
}
