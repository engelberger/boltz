# Mutation and Masking in Boltz

This document describes the mutation and masking functionality added to Boltz, which allows users to:

1. **Mutate** specific amino acids in protein sequences before prediction
2. **Mask** specific positions in the Multiple Sequence Alignment (MSA) during inference

## Features

### Mutations

Mutations allow you to change specific amino acids in the protein sequence. The mutated sequence will be used for:
- MSA generation (if using `--use_msa_server`)
- Structure prediction
- All downstream analysis

### Masking

Masking allows you to hide evolutionary information at specific positions in the MSA by:
- Replacing amino acids with a mask token (default: "X") 
- Optionally zeroing deletion matrix information
- Preserving the target sequence (first row of MSA is never masked)

## Command Line Usage

### Basic Mutation Examples

```bash
# Single mutation: Change Alanine at position 5 to Glycine
# Default: Only mutates target sequence, keeps original MSA for fair comparison
boltz predict input.fasta --mutations "A5G"

# Multiple mutations
boltz predict input.fasta --mutations "A5G,R10C,V20L"

# Chain-specific mutations (for multi-chain proteins)
boltz predict input.fasta --mutations "A:G10C,B:R20D"

# Regenerate MSA with mutated sequence (for studying evolutionary context)
boltz predict input.fasta --mutations "A5G" --mutate_msa_query
```

### Basic Masking Examples

```bash
# Mask specific positions in MSA
boltz predict input.fasta --mask_positions "10,15,20" --mask_msa

# Mask a range of positions
boltz predict input.fasta --mask_positions "10-15,20-25" --mask_msa

# Mask with deletion matrix zeroing
boltz predict input.fasta --mask_positions "10-15" --mask_msa --mask_deletion_matrix

# Use different mask token
boltz predict input.fasta --mask_positions "10-15" --mask_msa --mask_token "-"
```

### Combined Usage

```bash
# Apply mutations and masking together
boltz predict input.fasta \
  --mutations "A5G,R10C" \
  --mask_positions "15-20,25" \
  --mask_msa \
  --mask_deletion_matrix
```

## Command Line Options

### Mutation Options

- `--mutations TEXT`: Comma-separated list of mutations
  - Format: `[chain:]original_aa{position}new_aa`
  - Examples: `A123G`, `B:R45C`, `A5G,V10L`
  - Positions are 1-based (as typical in biology)

- `--mutate_msa_query`: Regenerate MSA using the mutated sequence
  - Default: Only mutate target sequence, keep original MSA
  - Use this flag to study evolutionary context of mutations

### Masking Options

- `--mask_positions TEXT`: Positions to mask in MSA
  - Format: Single positions (`10,20,30`) or ranges (`10-15,20-25`)
  - Positions are 1-based
  
- `--mask_msa`: Enable MSA masking at specified positions

- `--mask_deletion_matrix`: Enable deletion matrix masking at specified positions

- `--mask_token [A|R|N|D|C|E|Q|G|H|I|L|K|M|F|P|S|T|W|Y|V|X|-]`: Character to use for masking (default: X)

## Input Format Support

### FASTA Format

Mutations and masking are controlled entirely through command-line options when using FASTA input:

```fasta
>A|protein|path/to/msa.a3m
MVTPEGNVSLVDESLLVGVTDEDRAV
```

### YAML Format

For YAML input, mutations can be applied by directly modifying the sequence in the YAML file. Masking must still be controlled through command-line options:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVTPEGNVSLVDESLLVGVTDEDRAV  # Already mutated sequence
      msa: path/to/msa.a3m
```

## Technical Details

### Mutation Implementation

**Default Behavior (Target-Only Mutations):**
1. **Target Sequence Only**: By default, mutations are applied only to the query sequence (first sequence in MSA)
2. **MSA Preservation**: The rest of the MSA remains unchanged for fair wild-type vs mutant comparison
3. **Late Application**: Mutations are applied after MSA generation during tokenization

**Alternative Behavior (with `--mutate_msa_query`):**
1. **Early Application**: Mutations are applied immediately after parsing input, before MSA generation
2. **MSA Regeneration**: When using `--use_msa_server`, the mutated sequence is used for MSA search
3. **Full Context**: Studies how mutations affect evolutionary search and alignment

**Common Features:**
- **Multi-chain Support**: Chain-specific mutations use the format `chain:mutation`
- **Validation**: Original amino acid is checked and warnings are issued for mismatches

### Masking Implementation

1. **Late Application**: Masking is applied after tokenization but before featurization
2. **Target Preservation**: The first sequence (query) in the MSA is never masked
3. **Deletion Handling**: When `--mask_deletion_matrix` is used, deletion information at masked positions is zeroed
4. **Token Replacement**: Masked positions are replaced with the specified mask token (default: X/UNK)

### Position Indexing

- **User Input**: 1-based indexing (position 1 = first amino acid)
- **Internal**: 0-based indexing (position 0 = first amino acid)
- **Conversion**: Handled automatically by the parser

## Example Workflows

### Mutation Study

To study the effect of specific mutations (recommended approach for fair comparison):

```bash
# Wild type
boltz predict protein.fasta --out_dir wt_results

# Single mutant (target-only mutation, same MSA)
boltz predict protein.fasta --mutations "A123G" --out_dir mutant_A123G

# Double mutant (target-only mutation, same MSA)
boltz predict protein.fasta --mutations "A123G,R45C" --out_dir mutant_double

# Compare with regenerated MSA (different evolutionary context)
boltz predict protein.fasta --mutations "A123G" --mutate_msa_query --out_dir mutant_A123G_new_msa
```

### MSA Masking Experiment

To study the importance of specific positions:

```bash
# Full MSA
boltz predict protein.fasta --out_dir full_msa

# Mask active site region (positions 10-15)
boltz predict protein.fasta \
  --mask_positions "10-15" \
  --mask_msa \
  --out_dir masked_active_site
```

### Combined Analysis

Study mutant with reduced evolutionary information:

```bash
boltz predict protein.fasta \
  --mutations "A123G" \
  --mask_positions "120-130" \
  --mask_msa \
  --mask_deletion_matrix \
  --out_dir mutant_with_masking
```

## Error Handling

### Mutation Errors (Fatal)

- Invalid mutation format → Parse error with helpful message
- Position out of bounds → Error with sequence length information  
- Invalid amino acid codes → Error with valid amino acid list
- Chain not found → Error with available chain list

### Masking Errors (Non-fatal)

- Invalid positions → Warning, invalid positions skipped
- Mask token not recognized → Error with valid token list
- MSA loading failures → Warning, proceed without masking

## Implementation Notes

- **Default mutations**: Modify only the MSA query sequence (target-only) for fair comparison
- **With `--mutate_msa_query`**: Mutations modify the `Target` object before MSA generation
- **Masking**: Modifies the `Tokenized` object after MSA loading
- Both features are compatible with Boltz1 and Boltz2
- All changes are applied in-memory; original files are not modified
- Extensive logging provides feedback on applied changes

## Limitations

- YAML masking syntax is not yet implemented (use command-line options)
- Only amino acid substitutions are supported (no insertions/deletions)
- Chain-specific masking is not yet available
- Template masking is not supported

## Future Enhancements

- YAML syntax for masking configuration
- Chain-specific masking options
- Support for insertions and deletions
- Template structure masking
- Position-specific mask tokens 