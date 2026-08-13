"""Privacy-safe ingestion exceptions."""


class GenotypeParseError(ValueError):
    """A structural source-record failure without genotype content."""

    def __init__(self, line_number: int, reason: str) -> None:
        self.line_number = line_number
        self.reason = reason
        super().__init__(f"line {line_number}: {reason}")
