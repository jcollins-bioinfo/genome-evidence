"""Minimal DuckDB DDL for reconstructible analytical state."""

SCHEMA_VERSION = 1
DDL = """
CREATE TABLE IF NOT EXISTS schema_metadata(version INTEGER NOT NULL);
CREATE TABLE IF NOT EXISTS artifact_registry (
    artifact_id VARCHAR PRIMARY KEY,
    artifact_type VARCHAR NOT NULL,
    parquet_uri VARCHAR NOT NULL,
    checksum VARCHAR NOT NULL,
    run_id VARCHAR NOT NULL,
    created_at TIMESTAMP NOT NULL
);
"""
