"""Lexical-only genotype and chromosome categorization."""

from genome_evidence.qc.models import LexicalGenotypeCategory, LexicalZygosity

NO_CALL_TOKENS = {"--", "", "NC", "NOCALL"}
RECOGNIZED_CHROMOSOMES = {str(value) for value in range(1, 23)} | {"X", "Y", "MT", "M", "XY"}


def categorize_genotype(token: str) -> tuple[LexicalGenotypeCategory, LexicalZygosity | None]:
    if token.upper() in NO_CALL_TOKENS:
        return LexicalGenotypeCategory.NO_CALL, None
    if len(token) == 1:
        return LexicalGenotypeCategory.SINGLE_ALLELE_TOKEN, None
    if len(token) == 2:
        zygosity = None
        if set(token.upper()) <= set("ACGT"):
            zygosity = (
                LexicalZygosity.HOMOZYGOUS_LEXICAL
                if token[0].upper() == token[1].upper()
                else LexicalZygosity.HETEROZYGOUS_LEXICAL
            )
        return LexicalGenotypeCategory.TWO_ALLELE_TOKEN, zygosity
    return LexicalGenotypeCategory.OTHER_TOKEN, None
