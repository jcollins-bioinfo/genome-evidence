# ADR-0003: DuckDB and Parquet storage strategy

**Status:** Accepted

## Decision

Use DuckDB for local relational analytical/query state, Parquet for immutable or versioned tabular artifacts, and Pydantic models at Python API boundaries. Use Polars for dataframe-oriented transformations.

## Rationale and consequences

This combination is lightweight, reproducible, columnar, and appropriate for a personal local analytical system. PostgreSQL would add operational burden before concurrent service requirements exist. A graph database would duplicate the relational source of truth and encourage premature modelling; directed evidence relationships fit relational tables. Revisit those choices only when measured requirements justify them.
