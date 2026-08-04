# Literature Review Documentation

This directory contains analysis of published literature with associated RNA-seq and DNA methylation datasets for Crassostrea/Magallana gigas.

## Files in this directory:

- `key-publications.md` - Verified publications (real PMIDs + DOIs) with associated/related genomic datasets, organized by theme
- `pubmed-search-results.md` - Executed Phase 2 searches: queries, result counts, query URLs, journal/author distributions, and PubMed→archive dataset linkage (Phase 2 deliverable)

## Method

Phase 2 literature mining was executed live against the NCBI PubMed database
using the E-utilities API (`esearch`/`esummary`/`elink`), then cross-linked to
BioProject/SRA/GEO. Every publication in `key-publications.md` carries a verified
PMID and DOI. Dataset accessions are marked as verified only where confirmed
against NCBI; unverified availability is flagged as such.

## Focus Areas:
- Roberts Lab publications (University of Washington)
- International oyster genomics studies (France: Caen/Perpignan; China: IOCAS)
- Environmental / transgenerational epigenetics
- Developmental biology and reproduction
- Immunity and disease (POMS, Vibrio, OsHV-1)

## Time Period Covered: 2010-2026

## Last Updated: 2026-08-04 (Phase 2 executed)