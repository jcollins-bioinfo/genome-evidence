"""Versioned external-evidence ingestion and exact-linking public API."""

from .models import AnnotationConfig, ClinVarIngestionConfig
from .pipeline import ingest_clinvar_vcv, link_external_evidence

__all__ = [
    "AnnotationConfig",
    "ClinVarIngestionConfig",
    "ingest_clinvar_vcv",
    "link_external_evidence",
]
