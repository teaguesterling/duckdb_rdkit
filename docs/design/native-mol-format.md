# Design: Native RDKit Mol Format (v2)

## Summary

Introduce native RDKit molecule storage as the default `Mol` type, with a separate `DalkeFP` fingerprint type for optimized substructure screening. Keep `UmbraMol` as a legacy format for backwards compatibility.

## Type Architecture

### `Mol` (NEW - Default)
Pure RDKit MolPickler binary format.

```
┌─────────────────────────────────────────────────────┐
│         RDKit MolPickler Binary (variable)          │
└─────────────────────────────────────────────────────┘
```

- **Storage**: BLOB (PhysicalType::VARCHAR for pointer swizzling)
- **Compatibility**: Direct RDKit interoperability
- **Usage**: Default type for molecule storage

### `DalkeFP` (NEW)
64-bit substructure screening fingerprint (Dalke Extended).

```
┌────────────────────────────────────────────────────────────────────┐
│ Dalke Fragment Bits (55) │ HeavyAtoms(4) │ Rings(2) │ Flags(3)    │
├──────────────────────────┼───────────────┼──────────┼─────────────┤
│ Bits 0-54                │ Bits 55-58    │ Bits 59-60│ Bits 61-63 │
└────────────────────────────────────────────────────────────────────┘
```

**Bit Layout:**
| Bits | Field | Description |
|------|-------|-------------|
| 0-54 | Fragment bits | Original Dalke fragment/count patterns |
| 55-58 | Heavy atom bucket | 4 bits (16 size ranges) |
| 59-60 | Ring count | 2 bits (0, 1, 2, 3+) |
| 61 | Has stereocenters | 1 if molecule has chiral centers |
| 62 | Has charges | 1 if molecule has formal charges |
| 63 | Reserved | Future use |

**Heavy Atom Buckets (4 bits = 16 values):**
| Value | Range | Value | Range |
|-------|-------|-------|-------|
| 0 | 1-5 | 8 | 41-50 |
| 1 | 6-10 | 9 | 51-60 |
| 2 | 11-15 | 10 | 61-75 |
| 3 | 16-20 | 11 | 76-90 |
| 4 | 21-25 | 12 | 91-110 |
| 5 | 26-30 | 13 | 111-140 |
| 6 | 31-35 | 14 | 141-180 |
| 7 | 36-40 | 15 | 181+ |

**Substructure Screening Logic:**
```cpp
bool can_be_substruct(uint64_t target_fp, uint64_t query_fp) {
    // 1. Size check: target must be >= query size
    uint8_t target_size = (target_fp >> 55) & 0xF;
    uint8_t query_size = (query_fp >> 55) & 0xF;
    if (target_size < query_size) return false;

    // 2. Ring check: target must have >= query rings
    uint8_t target_rings = (target_fp >> 59) & 0x3;
    uint8_t query_rings = (query_fp >> 59) & 0x3;
    if (target_rings < query_rings) return false;

    // 3. Stereo check: if query has stereo, target must too
    if ((query_fp & (1ULL << 61)) && !(target_fp & (1ULL << 61))) return false;

    // 4. Charge check: if query has charges, target must too
    if ((query_fp & (1ULL << 62)) && !(target_fp & (1ULL << 62))) return false;

    // 5. Fragment bits: all query bits must be in target
    uint64_t frag_mask = (1ULL << 55) - 1;  // bits 0-54
    if ((target_fp & query_fp & frag_mask) != (query_fp & frag_mask)) return false;

    return true;  // Might be substruct, need full verification
}
```

- **Storage**: UINT64 (8 bytes, fixed)
- **Purpose**: Fast substructure search filtering with extended screening
- **Based on**: Andrew Dalke's fragment patterns + extended molecular properties

### `UmbraMol` (Single-Column Optimized)
Combined format with embedded fingerprint prefix.

```
┌─────────────────────────────────────────────────────┐
│ DalkeFP (8B) │ RDKit MolPickler Binary (variable)  │
└─────────────────────────────────────────────────────┘
```

- **Storage**: BLOB with 8-byte prefix
- **Purpose**: Single-column storage with embedded substructure filter
- **Use case**: When you want FP optimization without managing two columns
- **Trade-off**: Slightly larger scan (includes mol blob) vs separate FP column

---

## Functions

### Mol Functions
```sql
-- Construction
mol_from_smiles(VARCHAR) -> Mol
mol_from_ctab(VARCHAR) -> Mol          -- future
mol_from_pkl(BLOB) -> Mol              -- future

-- Conversion
mol_to_smiles(Mol) -> VARCHAR
mol_to_pkl(Mol) -> BLOB                -- just returns the blob
mol_to_umbramol(Mol) -> UmbraMol       -- convert to legacy format

-- Casts
'CCO'::Mol                             -- implicit mol_from_smiles
```

### DalkeFP Functions
```sql
-- Construction
dalke_fp(Mol) -> DalkeFP
dalke_fp(UmbraMol) -> DalkeFP          -- extract from UmbraMol

-- Operators (for substructure screening)
dalke_fp_contains(DalkeFP, DalkeFP) -> BOOLEAN  -- target contains query bits
fp1 @> fp2                                       -- operator alias

-- The contains check: (target & query) == query
-- All bits set in query must be set in target
```

### UmbraMol Functions
```sql
-- Construction
umbramol_from_smiles(VARCHAR) -> UmbraMol
umbramol_from_mol(Mol) -> UmbraMol

-- Extraction
umbramol_to_mol(UmbraMol) -> Mol
umbramol_get_fp(UmbraMol) -> DalkeFP

-- Casts
umbramol::Mol                          -- extract mol
umbramol::DalkeFP                      -- extract fingerprint
```

### Substructure Search
```sql
-- Works on all mol types
is_substruct(Mol, Mol) -> BOOLEAN
is_substruct(UmbraMol, Mol) -> BOOLEAN  -- uses embedded FP for speedup
is_substruct(UmbraMol, UmbraMol) -> BOOLEAN

-- Manual optimization with separate FP column
SELECT * FROM molecules
WHERE fp @> dalke_fp(query_mol)         -- fast FP filter
  AND is_substruct(mol, query_mol);     -- full verification
```

---

## Usage Patterns

### Pattern 1: Simple (New Default)
```sql
CREATE TABLE molecules (
    id INTEGER,
    mol MOL
);

INSERT INTO molecules
SELECT id, mol_from_smiles(smiles) FROM raw_data;

-- Substructure search (computes FP on-demand internally)
SELECT * FROM molecules
WHERE is_substruct(mol, 'c1ccccc1'::Mol);
```

### Pattern 2: Optimized with Separate FP Column
```sql
CREATE TABLE molecules (
    id INTEGER,
    mol MOL,
    fp DALKEFP
);

INSERT INTO molecules
SELECT id,
       mol_from_smiles(smiles) as mol,
       dalke_fp(mol_from_smiles(smiles)) as fp
FROM raw_data;

-- Fast substructure search (scans 8-byte FP column first)
SELECT * FROM molecules
WHERE fp @> dalke_fp('c1ccccc1'::Mol)
  AND is_substruct(mol, 'c1ccccc1'::Mol);
```

### Pattern 3: Legacy UmbraMol
```sql
CREATE TABLE molecules (
    id INTEGER,
    mol UMBRAMOL
);

INSERT INTO molecules
SELECT id, umbramol_from_smiles(smiles) FROM raw_data;

-- Automatic FP optimization (uses embedded FP)
SELECT * FROM molecules
WHERE is_substruct(mol, 'c1ccccc1'::Mol);
```

---

## Performance Characteristics

| Pattern | Storage/Row | Substructure Scan | Notes |
|---------|-------------|-------------------|-------|
| Mol only | ~N bytes | O(N) full scan | Simple, compatible |
| Mol + DalkeFP | ~N+8 bytes | O(8) FP scan + O(M) verify | Fastest for large tables |
| UmbraMol | ~N+8 bytes | O(N+8) scan with early exit | Legacy, single column |

Where:
- N = molecule pickle size (varies, typically 100-1000 bytes)
- M = number of FP-passing candidates

**Key insight**: Separate columns allow DuckDB to scan just the 8-byte FP column, skipping the large mol blobs entirely until verification.

---

## Implementation Plan

### Phase 1: Add New Types
1. [ ] Implement `DalkeFP` type (UINT64 wrapper)
2. [ ] Add `dalke_fp(Mol) -> DalkeFP` function
3. [ ] Add `dalke_fp_contains(DalkeFP, DalkeFP) -> BOOLEAN`

### Phase 2: Refactor Mol/UmbraMol
1. [ ] Rename current `Mol` type to `UmbraMol`
2. [ ] Create new `Mol` type (pure RDKit pickle)
3. [ ] Update `mol_from_smiles()` to return new `Mol`
4. [ ] Add `umbramol_from_smiles()` for UmbraMol creation

### Phase 3: Conversion Functions
1. [ ] Add `mol_to_umbramol(Mol) -> UmbraMol`
2. [ ] Add `umbramol_to_mol(UmbraMol) -> Mol`
3. [ ] Add `umbramol_to_struct(UmbraMol) -> STRUCT(mol Mol, fp DalkeFP)`
4. [ ] Add `dalke_fp(UmbraMol) -> DalkeFP` (extract from prefix)

### Phase 4: Update Functions
1. [ ] Update `is_substruct()` for both Mol and UmbraMol
2. [ ] Update `is_exact_match()` for both types
3. [ ] Update all descriptor functions for both types
4. [ ] Update `mol_to_smiles()` for both types

### Phase 5: Reader Support
1. [ ] Add `moltype` parameter to `read_sdf()`
2. [ ] Update SDF scanner to support both output types

---

## Future Considerations (Not in Scope)

### AnnotatedMol Type
A rich composite type for advanced use cases:
```sql
AnnotatedMol = STRUCT(
    mol MOL,
    struct_fp DALKEFP,
    metadata MAP(VARCHAR, VARCHAR),
    fps MAP(VARCHAR, BFP)
)
```

This would allow storing multiple fingerprint types and arbitrary metadata alongside molecules.

---

## Reader Functions

All reader functions support a `moltype` parameter:

```sql
-- SDF reading
read_sdf('file.sdf', COLUMNS={...})                    -- default: Mol
read_sdf('file.sdf', COLUMNS={...}, moltype := 'umbra') -- returns UmbraMol

-- SMILES reading (future)
read_smiles('file.smi', moltype := 'default')          -- returns Mol
read_smiles('file.smi', moltype := 'umbra')            -- returns UmbraMol

-- mol_from_smiles
mol_from_smiles('CCO')                                 -- returns Mol
umbramol_from_smiles('CCO')                            -- returns UmbraMol
```

---

## Conversion Functions

```sql
-- Mol <-> UmbraMol
mol_to_umbramol(Mol) -> UmbraMol           -- wrap mol with computed FP
umbramol_to_mol(UmbraMol) -> Mol           -- extract just the mol

-- UmbraMol decomposition
umbramol_to_struct(UmbraMol) -> STRUCT(mol Mol, fp DalkeFP)

-- Example usage
SELECT umbramol_to_struct(mol).mol,
       umbramol_to_struct(mol).fp
FROM umbra_table;
```

---

## Breaking Change

This is a **breaking change**. Existing databases with the old `Mol` type (which was UmbraMol internally) will need migration:

```sql
-- Old tables have UmbraMol data stored as "Mol"
-- After upgrade, that type becomes UmbraMol

-- To convert to new Mol format:
CREATE TABLE new_table AS
SELECT id, umbramol_to_mol(mol) as mol
FROM old_table;
```

**No auto-detection**: `Mol` is always pure RDKit pickle, `UmbraMol` is always prefixed format. Types are strict.
