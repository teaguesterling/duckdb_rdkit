# Issue: Additional File Format Scanners/Readers

**Priority: P4 (Medium)**
**Labels:** enhancement, feature

## Summary
Implement additional file format readers beyond SDF, focusing on commonly used molecular file formats.

## Proposed Formats

### SMILES Files (High Priority)
**Extensions:** `.smi`, `.smiles`, `.txt`

SMILES files are simple text files with one molecule per line:
```
CCO ethanol
c1ccccc1 benzene
CC(=O)O acetic_acid
```

Format variations:
- Tab-separated: `SMILES\tID\tprop1\tprop2`
- Space-separated: `SMILES ID`
- SMILES only (one per line)
- With header row or without

```sql
-- Proposed API
read_smiles('molecules.smi')
read_smiles('molecules.smi', COLUMNS={'smiles': 'MOL', 'name': 'VARCHAR'})
read_smiles('molecules.smi', delimiter := '\t', header := true)

-- Replacement scan
SELECT * FROM 'molecules.smi';
```

### MOL2 Files (Medium Priority)
**Extensions:** `.mol2`

Tripos MOL2 format, commonly used for:
- 3D coordinates
- Partial charges
- Protein-ligand docking

```sql
read_mol2('ligands.mol2')
read_mol2('ligands.mol2', COLUMNS={'mol': 'MOL', 'name': 'VARCHAR'})
```

### CSV with SMILES Column (Medium Priority)
Auto-detect SMILES columns in CSV files.

```sql
-- Explicit SMILES column
read_csv('data.csv', smiles_column := 'structure')

-- Auto-detect (look for columns named 'smiles', 'SMILES', 'structure', etc.)
read_csv_mol('data.csv')
```

### Compressed Files (Low Priority)
Support gzipped versions of all formats:
- `.sdf.gz`, `.smi.gz`, `.mol2.gz`

## Implementation Plan

### Phase 1: SMILES Reader
1. Create `src/smiles_scanner/` directory structure
2. Implement basic SMILES file parsing
3. Support tab and space delimiters
4. Support optional header row
5. Auto-detect column count from first few lines
6. Register replacement scan for `.smi` and `.smiles`

### Phase 2: MOL2 Reader
1. Use RDKit's Mol2MolSupplier
2. Handle multi-molecule MOL2 files
3. Extract common properties (name, charges)

### Phase 3: Enhanced CSV Support
1. Add SMILES column detection to CSV reader
2. Auto-convert detected columns to Mol type

## Example Implementation (SMILES)
```cpp
// Similar structure to SDF scanner
class SMILESScanner {
  void AutoDetect(...);  // Detect delimiter, header, columns
  void ExtractNextChunk(...);  // Parse lines into molecules
};
```

## Tasks
- [ ] Design SMILES scanner architecture
- [ ] Implement basic SMILES file reading
- [ ] Support delimiter options (tab, space, comma)
- [ ] Support header/no-header modes
- [ ] Add SMILES replacement scan
- [ ] Implement MOL2 scanner
- [ ] Add gzip support
- [ ] Write comprehensive tests
- [ ] Document all file format options
