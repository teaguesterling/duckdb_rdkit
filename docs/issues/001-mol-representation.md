# Issue: RDKit-native Mol Representation

**Priority: P1 (Highest)**
**Labels:** enhancement, architecture

## Summary
Evaluate and potentially migrate from the current "UmbraMol" format to a standard RDKit-native molecule representation.

## Current State
The current Mol type uses a custom "UmbraMol" format:
- **Prefix (8 bytes):** Dalke fingerprint for fast substructure filtering
  - 4 bytes inlined in DuckDB's string_t PREFIX
  - 4 bytes at the start of blob data
  - Total: 55-bit fingerprint
- **Body:** RDKit MolPickler binary format

## Issues with Current Format
1. **Proprietary format** - Cannot be directly read by other RDKit tools
2. **Requires conversion** - `mol_to_rdkit_mol()` needed for interop
3. **Complexity** - Pointer swizzling for disk/memory transitions
4. **Maintenance burden** - Custom serialization code

## Benefits of Current Format
1. **Fast substructure filtering** - Dalke fingerprint enables early bailout
2. **Compact storage** - Fingerprint is only 8 bytes overhead

## Proposed Options

### Option A: Keep UmbraMol (Status Quo)
- Pro: Fast substructure filtering
- Pro: No migration needed
- Con: Proprietary format, maintenance burden

### Option B: Standard RDKit Pickle
- Pro: Direct RDKit compatibility
- Pro: Simpler codebase
- Con: Lose fast filtering (compute fingerprint on-demand)
- Con: Breaking change for existing databases

### Option C: Hybrid Approach
- Store standard RDKit pickle as primary format
- Compute/cache fingerprint lazily or at query time
- Pro: Compatibility + can optimize later
- Con: More complex than B, may lose some performance

## Decision Criteria
- Performance benchmarks on substructure search
- Storage overhead comparison
- User feedback on interoperability needs
- Complexity of maintaining UmbraMol

## Tasks
- [ ] Benchmark substructure search with/without fingerprint prefix
- [ ] Measure storage overhead
- [ ] Prototype Option B and C
- [ ] Document migration path if format changes
- [ ] Make decision based on data

## Related
- Affects: All Mol storage and operations
- Blocks: Some fingerprint features (if we want consistent storage)
