class InvalidDate(Exception):
    """Custom exception class to handle dates in unsupported formats."""

    def __init__(self, date, message):
        """Initialize the error with the culprit date string and a meaningful message."""
        self.date = date
        self.message = message

    def __str__(self):
        """Return a human-readable summary of the error."""
        return f"Invalid date {self.date!r}: {self.message}"

class InvalidYearBounds(Exception):
    """Custom exception class to handle year bounds in unsupported formats."""
