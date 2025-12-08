# Issue: PostgreSQL RDKit Extension Feature Parity

**Priority: P3 (Medium-High)**
**Labels:** enhancement, compatibility

## Summary
Implement missing functions to achieve feature parity with the PostgreSQL RDKit cartridge, enabling users to migrate existing workflows.

## Current Status

### Implemented (10 functions)
- `mol_from_smiles()`, `mol_to_smiles()`
- `mol_to_rdkit_mol()`
- `is_substruct()`, `is_exact_match()`
- `mol_amw()`, `mol_exactmw()`, `mol_logp()`, `mol_tpsa()`
- `mol_hba()`, `mol_hbd()`, `mol_num_rotatable_bonds()`
- `mol_qed()`

### Missing - High Priority

#### Format Conversion
```sql
mol_from_ctab(text) -> mol          -- V2000/V3000 mol block
mol_to_ctab(mol) -> text
mol_from_pkl(bytea) -> mol          -- RDKit binary pickle
mol_to_pkl(mol) -> bytea
mol_to_smarts(mol) -> text          -- SMARTS pattern
mol_to_svg(mol, width, height) -> text  -- SVG visualization
mol_from_json(text) -> mol          -- RDKit JSON
mol_to_json(mol) -> text
```

#### Validation Functions
```sql
is_valid_smiles(text) -> boolean
is_valid_smarts(text) -> boolean
is_valid_ctab(text) -> boolean
```

#### Search Functions
```sql
substruct_count(mol, pattern) -> integer  -- Count matches
```

### Missing - Medium Priority

#### Integer Descriptors
```sql
mol_numatoms(mol) -> integer
mol_numheavyatoms(mol) -> integer
mol_numheteroatoms(mol) -> integer
mol_numrings(mol) -> integer
mol_numaromaticrings(mol) -> integer
mol_numaliphaticrings(mol) -> integer
mol_numsaturatedrings(mol) -> integer
mol_numaromaticheterocycles(mol) -> integer
mol_numaliphaticheterocycles(mol) -> integer
mol_numsaturatedheterocycles(mol) -> integer
mol_numspiroatoms(mol) -> integer
mol_numbridgeheadatoms(mol) -> integer
mol_numamidebonds(mol) -> integer
```

#### Float Descriptors
```sql
mol_fractioncsp3(mol) -> float
mol_labuteasa(mol) -> float
mol_chi0n(mol) -> float
mol_chi1n(mol) -> float
mol_chi2n(mol) -> float
mol_chi3n(mol) -> float
mol_chi4n(mol) -> float
mol_chi0v(mol) -> float
mol_chi1v(mol) -> float
mol_chi2v(mol) -> float
mol_chi3v(mol) -> float
mol_chi4v(mol) -> float
mol_kappa1(mol) -> float
mol_kappa2(mol) -> float
mol_kappa3(mol) -> float
mol_phi(mol) -> float
mol_hallkieralpha(mol) -> float
```

#### Text Descriptors
```sql
mol_formula(mol) -> text
mol_inchi(mol) -> text              -- Requires InChI support
mol_inchikey(mol) -> text           -- Requires InChI support
mol_murckoscaffold(mol) -> text
```

### Missing - Lower Priority

#### Query Molecule Type
```sql
-- New type for SMARTS-based queries
qmol_from_smarts(text) -> qmol
qmol_from_ctab(text) -> qmol
```

#### Reaction Support
```sql
-- New reaction type and functions
reaction_from_smiles(text) -> reaction
reaction_from_smarts(text) -> reaction
reaction_to_smiles(reaction) -> text
-- etc.
```

## Implementation Strategy
1. **Phase 1:** Format conversion functions (enables interop)
2. **Phase 2:** Additional descriptors (quick wins, use existing RDKit)
3. **Phase 3:** Query molecule type (enables advanced searching)
4. **Phase 4:** Reaction support (new functionality)

## Tasks
- [ ] Implement mol_from_ctab/mol_to_ctab
- [ ] Implement mol_to_svg
- [ ] Implement validation functions
- [ ] Implement all integer descriptors
- [ ] Implement all float descriptors
- [ ] Implement mol_formula
- [ ] Consider InChI support (requires additional RDKit build option)
- [ ] Add comprehensive tests
- [ ] Document PostgreSQL migration guide
