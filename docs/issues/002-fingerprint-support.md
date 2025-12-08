# Issue: Fingerprint Support (bfp/sfp types and similarity functions)

**Priority: P2 (High)**
**Labels:** enhancement, feature

## Summary
Implement fingerprint types and similarity functions, which are essential for cheminformatics workflows like similarity searching, clustering, and virtual screening.

## Required Types

### `bfp` - Bit Fingerprint
Fixed-size bit vector fingerprint.
- Storage: BLOB with fixed size header
- Operations: AND, OR, XOR, popcount, bit access

### `sfp` - Sparse Fingerprint
Count-based dictionary fingerprint (feature -> count).
- Storage: Variable-size, sorted (feature_id, count) pairs
- Operations: Add, subtract, intersection, union

## Required Fingerprint Generators

### Bit Fingerprints (-> bfp)
```sql
-- Morgan/Circular fingerprints (most widely used)
morganbv_fp(mol, radius := 2, nbits := 2048) -> bfp
featmorganbv_fp(mol, radius := 2, nbits := 2048) -> bfp

-- RDKit topological fingerprint
rdkit_fp(mol, minPath := 1, maxPath := 7, nbits := 2048) -> bfp

-- MACCS keys (166 predefined structural keys)
maccs_fp(mol) -> bfp

-- Atom pair and torsion
atompairbv_fp(mol, nbits := 2048) -> bfp
torsionbv_fp(mol, nbits := 2048) -> bfp
```

### Sparse Fingerprints (-> sfp)
```sql
-- Morgan with counts
morgan_fp(mol, radius := 2) -> sfp
featmorgan_fp(mol, radius := 2) -> sfp

-- Atom pair and torsion with counts
atompair_fp(mol) -> sfp
torsion_fp(mol) -> sfp
```

## Required Similarity Functions
```sql
-- Tanimoto/Jaccard coefficient (most common)
tanimoto_sml(fp1, fp2) -> FLOAT  -- Works with both bfp and sfp

-- Dice coefficient
dice_sml(fp1, fp2) -> FLOAT

-- Tversky index (asymmetric similarity)
tversky_sml(fp1, fp2, alpha, beta) -> FLOAT
```

## Required Utility Functions
```sql
-- Fingerprint size/statistics
bfp_size(bfp) -> INTEGER          -- Number of bits
bfp_popcount(bfp) -> INTEGER      -- Number of set bits
sfp_size(sfp) -> INTEGER          -- Number of features

-- Fingerprint manipulation
bfp_and(bfp, bfp) -> bfp
bfp_or(bfp, bfp) -> bfp
sfp_add(sfp, sfp) -> sfp
sfp_subtract(sfp, sfp) -> sfp
```

## Example Usage
```sql
-- Find similar molecules
SELECT mol_to_smiles(mol), tanimoto_sml(morganbv_fp(mol), morganbv_fp('CCO'::mol)) as sim
FROM molecules
WHERE tanimoto_sml(morganbv_fp(mol), morganbv_fp('CCO'::mol)) > 0.7
ORDER BY sim DESC;

-- Cluster molecules by fingerprint
SELECT mol_to_smiles(mol), maccs_fp(mol) as fp
FROM molecules;
```

## Implementation Notes
1. Consider using DuckDB's native BLOB type for storage
2. Implement efficient comparison operators for index support
3. May need custom aggregate functions for clustering
4. Consider lazy computation with caching for frequently used fingerprints

## PostgreSQL RDKit Reference
- https://www.rdkit.org/docs/Cartridge.html#fingerprint-functions
- https://www.rdkit.org/docs/Cartridge.html#similarity-functions

## Tasks
- [ ] Design bfp storage format
- [ ] Design sfp storage format
- [ ] Implement Morgan fingerprint generators
- [ ] Implement MACCS fingerprint
- [ ] Implement RDKit fingerprint
- [ ] Implement Tanimoto similarity
- [ ] Implement Dice similarity
- [ ] Add tests for all functions
- [ ] Add documentation
