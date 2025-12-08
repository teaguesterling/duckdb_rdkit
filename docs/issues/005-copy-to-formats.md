# Issue: COPY TO / Export Format Support

**Priority: P5 (Medium-Low)**
**Labels:** enhancement, feature

## Summary
Implement the ability to export molecular data to common file formats using DuckDB's COPY TO syntax or dedicated functions.

## Proposed Formats

### SDF Export (High Priority)
```sql
-- Export to SDF file
COPY molecules TO 'output.sdf' (FORMAT sdf);

-- With specific columns as properties
COPY (SELECT mol, name, mw FROM molecules) TO 'output.sdf' (FORMAT sdf);

-- Function-based export
SELECT write_sdf(mol, name, mw) FROM molecules;  -- Returns SDF text
```

### SMILES Export (High Priority)
```sql
-- Export to SMILES file
COPY molecules TO 'output.smi' (FORMAT smiles);

-- With custom columns
COPY (SELECT mol, id FROM molecules) TO 'output.smi' (FORMAT smiles);

-- Tab-separated with header
COPY molecules TO 'output.smi' (FORMAT smiles, HEADER true, DELIMITER '\t');
```

### MOL2 Export (Medium Priority)
```sql
COPY molecules TO 'output.mol2' (FORMAT mol2);
```

### CSV with SMILES (Medium Priority)
```sql
-- Standard CSV but mol columns become SMILES strings
COPY molecules TO 'output.csv' (FORMAT csv);
```

## Implementation Approach

### Option A: COPY TO Extension
Extend DuckDB's COPY TO with new format handlers.

Pros:
- Native DuckDB syntax
- Integrates with existing workflows

Cons:
- More complex to implement
- May require DuckDB internals knowledge

### Option B: Table Functions
Create functions that return formatted text.

```sql
-- Write to file via DuckDB's file system
COPY (SELECT sdf_record(mol, {'name': name, 'mw': mw}) FROM molecules)
TO 'output.sdf' (FORMAT csv, HEADER false, QUOTE '');
```

Pros:
- Simpler implementation
- More flexible

Cons:
- Less elegant syntax
- May have escaping issues

### Option C: Dedicated Export Functions
```sql
-- Returns nothing, writes directly to file
SELECT export_sdf('output.sdf', mol, name, mw) FROM molecules;

-- Or aggregate function
SELECT write_sdf_file('output.sdf', mol, struct_pack(name, mw)) FROM molecules;
```

## SDF Format Details
```
<molecule name>
  <software line>

<atom block>
<bond block>
M  END
> <property1>
value1

> <property2>
value2

$$$$
```

## Tasks
- [ ] Decide on implementation approach (A, B, or C)
- [ ] Implement SDF export
- [ ] Implement SMILES export
- [ ] Handle Mol -> SMILES conversion for CSV
- [ ] Add compression support (.gz)
- [ ] Add MOL2 export
- [ ] Write tests
- [ ] Document export options

## Dependencies
- May benefit from having standard Mol format (Issue #001)
