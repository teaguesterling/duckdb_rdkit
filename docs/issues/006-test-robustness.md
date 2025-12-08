# Issue: Test Suite Robustness

**Priority: P6 (Ongoing)**
**Labels:** testing, quality

## Summary
Improve test coverage, edge case handling, and overall quality assurance of the extension.

## Current Test Coverage

### Existing Test Files
- `test/sql/mol_conversion.test` - Mol type conversion
- `test/sql/mol_search.test` - Substructure and exact match
- `test/sql/mol_descriptors.test` - Descriptor functions
- `test/sql/mol_storage.test` - Storage and retrieval
- `test/sql/edge_cases.test` - Edge cases
- `test/sql/qed.test` - QED calculation
- `test/sql/sdf_scanner/read_sdf.test` - SDF reading
- `test/sql/sdf_scanner/sdf_replacement_scan.test` - SDF auto-detection

### Current Stats
- 8 test files
- ~242 assertions (as of latest)

## Areas for Improvement

### 1. Edge Case Testing
- [ ] Empty inputs (empty strings, NULL values)
- [ ] Invalid SMILES strings (various malformed inputs)
- [ ] Very large molecules (performance/memory)
- [ ] Unicode in molecule names/properties
- [ ] Special characters in file paths
- [ ] Empty files (done for SDF)
- [ ] Malformed SDF files (missing delimiters, truncated)

### 2. Error Message Testing
- [ ] Verify error messages are helpful and actionable
- [ ] Test error conditions don't cause crashes
- [ ] Ensure proper cleanup on error

### 3. Performance Testing
- [ ] Benchmark substructure search on large datasets
- [ ] Test memory usage with large molecules
- [ ] Profile common operations
- [ ] Add performance regression tests

### 4. Cross-Platform Testing
- [ ] Linux (primary, CI)
- [ ] macOS (CI)
- [ ] Windows (needs attention)

### 5. Integration Testing
- [ ] Multiple operations in same query
- [ ] Joins between mol tables
- [ ] Aggregations with mol columns
- [ ] Persistence (save/load database)

### 6. Concurrency Testing
- [ ] Parallel queries
- [ ] Thread safety of RDKit operations

## Specific Test Cases to Add

### SMILES Parsing Edge Cases
```sql
-- Should handle or error gracefully:
SELECT mol_from_smiles('');           -- Empty string
SELECT mol_from_smiles('   ');        -- Whitespace only
SELECT mol_from_smiles('invalid');    -- Invalid SMILES
SELECT mol_from_smiles('C' || repeat('C', 10000));  -- Very long
SELECT mol_from_smiles('[invalid]');  -- Invalid atom
```

### Descriptor Edge Cases
```sql
-- Single atom molecules
SELECT mol_logp('C'::mol);
SELECT mol_tpsa('[Na+]'::mol);

-- Disconnected fragments
SELECT mol_amw('CCO.CC'::mol);

-- Charged molecules
SELECT mol_hba('[NH4+]'::mol);
```

### SDF Scanner Edge Cases
```sql
-- Various malformed files
FROM read_sdf('truncated.sdf');
FROM read_sdf('missing_delimiter.sdf');
FROM read_sdf('binary_garbage.sdf');
```

## Testing Infrastructure

### Improvements Needed
- [ ] Add test data generation scripts
- [ ] Create performance benchmark suite
- [ ] Add fuzzing for SMILES parsing
- [ ] Set up continuous benchmarking

### Test Data
- [ ] Create canonical test SDF files
- [ ] Add large molecule test cases
- [ ] Include known edge case molecules

## Tasks
- [ ] Audit current test coverage
- [ ] Add missing edge case tests
- [ ] Improve error message tests
- [ ] Add performance benchmarks
- [ ] Set up test data management
- [ ] Document testing best practices
