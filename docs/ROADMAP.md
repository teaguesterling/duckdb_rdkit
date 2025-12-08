# DuckDB RDKit Extension Roadmap

## Prioritization Analysis

This document outlines the strategic priorities for the duckdb_rdkit extension based on:
- Foundational impact (changes that affect other features)
- User value (common use cases in cheminformatics)
- Implementation complexity
- PostgreSQL RDKit extension parity

## Priority Order

### Priority 1: RDKit-native Mol Representation
**Impact: Foundational | Effort: Medium-High**

Currently, we use a custom "UmbraMol" format that embeds a Dalke fingerprint prefix for fast substructure filtering. While this provides performance benefits, it creates a proprietary format that:
- Cannot be directly read by other RDKit tools
- Requires conversion functions (`mol_to_rdkit_mol`)
- May complicate future feature development

**Decision needed:** Whether to:
1. Keep UmbraMol for performance benefits
2. Switch to standard RDKit pickle format
3. Hybrid approach (store both, or compute fingerprint on-demand)

### Priority 2: Fingerprint Support
**Impact: High User Value | Effort: High**

Fingerprints are essential for:
- Similarity searching (find similar molecules)
- Clustering and diversity analysis
- Virtual screening
- Machine learning feature generation

**Required types:**
- `bfp` - Bit fingerprint (fixed-size bit vector)
- `sfp` - Sparse fingerprint (count-based)

**Required fingerprint generators:**
- Morgan/Circular fingerprints (most common)
- RDKit fingerprint
- MACCS keys (166 bits, interpretable)
- Atom pair fingerprints
- Topological torsion fingerprints

**Required similarity functions:**
- `tanimoto_sml(fp, fp)` - Tanimoto/Jaccard similarity
- `dice_sml(fp, fp)` - Dice coefficient
- `tversky_sml(fp, fp, a, b)` - Tversky index

### Priority 3: PostgreSQL Extension Feature Parity
**Impact: Migration Path | Effort: High**

Users migrating from PostgreSQL RDKit need compatible functions. Current gaps:

**Missing Descriptors (25+):**
- `mol_formula()`, `mol_inchi()`, `mol_inchikey()`
- `mol_numheavyatoms()`, `mol_numheteroatoms()`, `mol_numrings()`
- `mol_numaromaticrings()`, `mol_numaliphaticrings()`
- `mol_fractioncsp3()`, `mol_labuteasa()`
- Chi connectivity indices (`mol_chi0n()` through `mol_chi4n()`)
- Kappa indices (`mol_kappa1()`, `mol_kappa2()`, `mol_kappa3()`)

**Missing Format Functions:**
- `mol_from_ctab()`, `mol_to_ctab()` - V2000/V3000 mol blocks
- `mol_from_json()`, `mol_to_json()` - RDKit JSON format
- `mol_to_smarts()` - SMARTS output
- `mol_to_svg()` - SVG visualization
- `is_valid_smiles()`, `is_valid_smarts()` - Validation

**Missing Search Functions:**
- `substruct_count()` - Count matches
- `substruct_chiral()` - Chirality-aware search

### Priority 4: Additional Scanners/Readers
**Impact: Data Ingestion | Effort: Medium**

Support reading common molecular file formats:
- **SMILES files** (.smi, .smiles) - Tab/space-separated SMILES + ID
- **MOL2 files** (.mol2) - Tripos format with 3D coordinates
- **CSV with SMILES column** - Auto-detect SMILES columns

### Priority 5: Copy-to/Export Formats
**Impact: Data Export | Effort: Medium**

Support writing molecular data to files:
- `COPY table TO 'file.sdf'` - SDF export
- `COPY table TO 'file.smi'` - SMILES export
- `COPY table TO 'file.mol2'` - MOL2 export

### Priority 6: Test Suite Robustness
**Impact: Quality | Effort: Ongoing**

Improve test coverage for:
- Edge cases (empty files, malformed input)
- Large datasets (performance testing)
- All supported platforms
- Regression tests for bug fixes

---

## Current Implementation Status

### Implemented Functions
| Category | Function | Notes |
|----------|----------|-------|
| Conversion | `mol_from_smiles()` | SMILES to Mol |
| Conversion | `mol_to_smiles()` | Mol to SMILES |
| Conversion | `mol_to_rdkit_mol()` | Mol to RDKit binary |
| Search | `is_exact_match()` | Exact structure match |
| Search | `is_substruct()` | Substructure search |
| Descriptor | `mol_amw()` | Molecular weight (approx) |
| Descriptor | `mol_exactmw()` | Molecular weight (exact) |
| Descriptor | `mol_logp()` | Wildman-Crippen LogP |
| Descriptor | `mol_tpsa()` | Topological PSA |
| Descriptor | `mol_qed()` | Drug-likeness |
| Descriptor | `mol_hba()` | H-bond acceptors |
| Descriptor | `mol_hbd()` | H-bond donors |
| Descriptor | `mol_num_rotatable_bonds()` | Rotatable bonds |
| I/O | `read_sdf()` | SDF file reader |
| I/O | `read_sdf_auto()` | SDF with auto-detect |

### PostgreSQL Parity Checklist

#### Types
- [x] `mol` - Molecule type
- [ ] `qmol` - Query molecule type
- [ ] `bfp` - Bit fingerprint type
- [ ] `sfp` - Sparse fingerprint type
- [ ] `reaction` - Chemical reaction type

#### Core Functions (High Priority)
- [x] `mol_from_smiles()`
- [x] `mol_to_smiles()`
- [ ] `mol_from_ctab()`
- [ ] `mol_to_ctab()`
- [ ] `mol_from_pkl()` / `mol_to_pkl()`
- [ ] `mol_to_svg()`
- [x] `substruct()` (as `is_substruct`)
- [ ] `substruct_count()`

#### Fingerprints (High Priority)
- [ ] `morgan_fp()` / `morganbv_fp()`
- [ ] `rdkit_fp()`
- [ ] `maccs_fp()`
- [ ] `tanimoto_sml()`
- [ ] `dice_sml()`

#### Descriptors
- [x] `mol_amw()`
- [x] `mol_exactmw()`
- [x] `mol_logp()`
- [x] `mol_tpsa()`
- [x] `mol_hba()`
- [x] `mol_hbd()`
- [x] `mol_numrotatablebonds()`
- [ ] `mol_formula()`
- [ ] `mol_inchi()`
- [ ] `mol_inchikey()`
- [ ] `mol_numheavyatoms()`
- [ ] `mol_numrings()`
- [ ] `mol_numaromaticrings()`
- [ ] ... (20+ more)
