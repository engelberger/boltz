"""Mutation and masking utilities for Boltz."""

import logging
import re
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from boltz.data import const
from boltz.data.types import MSA, Target, Tokenized


def parse_mutation_string(mutation_str: str) -> List[Tuple[Optional[str], int, str, str]]:
    """Parse mutation strings like 'A123G' or 'B:R45C' into structured data.
    
    Parameters
    ----------
    mutation_str : str
        Comma-separated mutation string (e.g., 'A123G,R45C' or 'A:G10C,B:R20D')
        
    Returns
    -------
    List[Tuple[Optional[str], int, str, str]]
        List of tuples: (chain_id, zero_based_position, original_aa, new_aa)
        chain_id is None if not specified
    """
    if not mutation_str or mutation_str.strip() == "":
        return []
        
    mutations = []
    for mut in mutation_str.split(","):
        mut = mut.strip()
        if not mut:
            continue
            
        # Pattern: optional chain_id followed by colon, then original_aa, position, new_aa
        # Examples: A123G, B:R45C, :A123G
        pattern = r"^(?:([A-Za-z0-9_]*):)?([A-Z])(\d+)([A-Z])$"
        match = re.match(pattern, mut)
        
        if not match:
            raise ValueError(f"Invalid mutation format: {mut}. Expected format: [chain:]original_aa{position}new_aa (e.g., A123G or B:R45C)")
            
        chain_id, original_aa, position_str, new_aa = match.groups()
        
        # Convert to 0-based indexing
        position = int(position_str) - 1
        if position < 0:
            raise ValueError(f"Position must be >= 1, got {position_str} in mutation {mut}")
            
        # Validate amino acids
        if original_aa not in const.prot_letter_to_token:
            raise ValueError(f"Invalid original amino acid: {original_aa} in mutation {mut}")
        if new_aa not in const.prot_letter_to_token:
            raise ValueError(f"Invalid new amino acid: {new_aa} in mutation {mut}")
            
        mutations.append((chain_id, position, original_aa, new_aa))
        
    return mutations


def parse_masking_positions(mask_str: str) -> List[int]:
    """Parse position strings like '10-15,20,30' into 0-based indices.
    
    Parameters
    ----------
    mask_str : str
        Position string with ranges and individual positions
        
    Returns
    -------
    List[int]
        List of 0-based indices to mask
    """
    if not mask_str or mask_str.strip() == "":
        return []
        
    positions = set()
    for part in mask_str.split(","):
        part = part.strip()
        if not part:
            continue
            
        if "-" in part:
            # Range like "10-15"
            try:
                start_str, end_str = part.split("-", 1)
                start = int(start_str.strip()) - 1  # Convert to 0-based
                end = int(end_str.strip()) - 1      # Convert to 0-based
                if start < 0 or end < 0:
                    raise ValueError(f"Positions must be >= 1, got range {part}")
                if start > end:
                    raise ValueError(f"Invalid range {part}: start > end")
                positions.update(range(start, end + 1))
            except ValueError as e:
                if "invalid literal" in str(e):
                    raise ValueError(f"Invalid range format: {part}. Expected format: start-end")
                raise
        else:
            # Single position like "20"
            try:
                pos = int(part) - 1  # Convert to 0-based
                if pos < 0:
                    raise ValueError(f"Position must be >= 1, got {part}")
                positions.add(pos)
            except ValueError:
                raise ValueError(f"Invalid position: {part}. Expected integer")
                
    return sorted(list(positions))


def apply_mutations_to_target(
    target: Target, 
    mutations: List[Tuple[Optional[str], int, str, str]]
) -> Target:
    """Apply mutations to the target sequences before MSA generation.
    
    Parameters
    ----------
    target : Target
        The target object with sequences to mutate
    mutations : List[Tuple[Optional[str], int, str, str]]
        List of mutations: (chain_id, position, original_aa, new_aa)
        
    Returns
    -------
    Target
        Target with mutated sequences
    """
    if not mutations:
        return target
        
    # Create a copy of sequences dict
    if target.sequences is None:
        logging.warning("No sequences found in target, skipping mutations")
        return target
        
    mutated_sequences = dict(target.sequences)
    
    # Group mutations by chain
    chain_mutations = {}
    for chain_id, position, original_aa, new_aa in mutations:
        if chain_id not in chain_mutations:
            chain_mutations[chain_id] = []
        chain_mutations[chain_id].append((position, original_aa, new_aa))
        
    # Apply mutations to each chain
    for chain_id, chain_mutations_list in chain_mutations.items():
        if chain_id is None:
            # Apply to all protein chains if no chain specified
            for entity_id, sequence in mutated_sequences.items():
                mutated_sequences[entity_id] = _apply_mutations_to_sequence(
                    sequence, chain_mutations_list, f"entity_{entity_id}"
                )
        else:
            # Find the entity_id for this chain_id
            entity_id = None
            for record_chain in target.record.chains:
                if record_chain.chain_name == chain_id:
                    entity_id = record_chain.entity_id
                    break
                    
            if entity_id is None:
                raise ValueError(f"Chain {chain_id} not found in target")
                
            if entity_id not in mutated_sequences:
                raise ValueError(f"No sequence found for chain {chain_id} (entity {entity_id})")
                
            mutated_sequences[entity_id] = _apply_mutations_to_sequence(
                mutated_sequences[entity_id], chain_mutations_list, chain_id
            )
    
    # Create new target with mutated sequences
    new_target = Target(
        record=target.record,
        structure=target.structure,
        sequences=mutated_sequences,
        residue_constraints=target.residue_constraints,
        templates=target.templates,
        extra_mols=target.extra_mols,
    )
    
    return new_target


def _apply_mutations_to_sequence(
    sequence: str, 
    mutations: List[Tuple[int, str, str]], 
    chain_name: str
) -> str:
    """Apply mutations to a single sequence string."""
    sequence_list = list(sequence)
    
    for position, original_aa, new_aa in mutations:
        if position >= len(sequence):
            raise ValueError(f"Position {position + 1} is out of bounds for chain {chain_name} (length {len(sequence)})")
            
        current_aa = sequence_list[position]
        if current_aa != original_aa:
            logging.warning(
                f"Chain {chain_name} position {position + 1}: expected {original_aa}, "
                f"found {current_aa}. Applying mutation anyway."
            )
            
        sequence_list[position] = new_aa
        logging.info(f"Applied mutation {chain_name}:{current_aa}{position + 1}{new_aa}")
        
    return "".join(sequence_list)


def apply_mutations_to_tokenized(
    tokenized: Tokenized,
    mutations: List[Tuple[Optional[str], int, str, str]]
) -> Tokenized:
    """Apply mutations to MSA query sequence only (first sequence in MSA).
    
    This function applies mutations only to the query sequence (target) 
    while keeping the rest of the MSA unchanged for fair comparison.
    
    Parameters
    ----------
    tokenized : Tokenized
        The tokenized structure with MSA
    mutations : List[Tuple[Optional[str], int, str, str]]
        List of mutations: (chain_id, position, original_aa, new_aa)
        
    Returns
    -------
    Tokenized
        Tokenized structure with mutations applied to MSA query only
    """
    if not mutations:
        return tokenized
        
    # If no MSA, fall back to token-level mutations
    if tokenized.msa is None:
        return _apply_mutations_to_tokens_fallback(tokenized, mutations)
    
    # Create a deep copy to avoid modifying original
    import copy
    result = copy.deepcopy(tokenized)
    
    # Apply mutations to MSA query sequence (first sequence only)
    for msa_chain_id, msa in result.msa.items():
        if len(msa.sequences) == 0:
            continue
            
        # Get the query sequence (first sequence)
        query_seq = msa.sequences[0]
        res_start, res_end = query_seq["res_start"], query_seq["res_end"]
        
        # Create mutable copy of residues
        new_residues = msa.residues.copy()
        
        # Apply mutations to query sequence only
        for chain_id, position, original_aa, new_aa in mutations:
            # Check if this mutation applies to this MSA chain
            # For now, apply to all chains if chain_id is None, or if chain matches
            if chain_id is not None and chain_id != msa_chain_id:
                continue
                
            # Convert position to 0-based and validate
            pos_0based = position - 1
            if pos_0based < 0 or pos_0based >= (res_end - res_start):
                logging.warning(f"Position {position} out of range for chain {msa_chain_id}")
                continue
                
            # Get new token ID
            if new_aa not in const.prot_letter_to_token:
                logging.warning(f"Unknown amino acid {new_aa}, skipping")
                continue
                
            new_token = const.prot_letter_to_token[new_aa]
            new_token_id = const.token_ids[new_token]
            
            # Apply mutation to query sequence
            global_pos = res_start + pos_0based
            if global_pos < len(new_residues):
                old_token_id = new_residues[global_pos]
                new_residues[global_pos] = new_token_id
                
                # Log the mutation
                old_aa = None
                for letter, token in const.prot_letter_to_token.items():
                    if const.token_ids[token] == old_token_id:
                        old_aa = letter
                        break
                        
                logging.info(f"Applied MSA query mutation: {old_aa or '?'}{position}{new_aa} in chain {msa_chain_id}")
        
        # Update MSA with new residues
        result.msa[msa_chain_id] = MSA(
            residues=new_residues,
            deletions=msa.deletions,
            sequences=msa.sequences,
        )
    
    return result


def _apply_mutations_to_tokens_fallback(
    tokenized: Tokenized,
    mutations: List[Tuple[Optional[str], int, str, str]]
) -> Tokenized:
    """Fallback: Apply mutations to tokenized structure (for when no MSA available).
    
    Parameters
    ----------
    tokenized : Tokenized
        The tokenized structure
    mutations : List[Tuple[Optional[str], int, str, str]]
        List of mutations: (chain_id, position, original_aa, new_aa)
        
    Returns
    -------
    Tokenized
        Tokenized structure with mutations applied
    """
    # Create mutable copy of tokens
    tokens = tokenized.tokens.copy()
    
    # Group mutations by chain
    chain_mutations = {}
    for chain_id, position, original_aa, new_aa in mutations:
        if chain_id not in chain_mutations:
            chain_mutations[chain_id] = []
        chain_mutations[chain_id].append((position, original_aa, new_aa))
    
    # Apply mutations
    for chain_id, chain_mutations_list in chain_mutations.items():
        if chain_id is None:
            # Apply to all protein chains
            unique_chains = np.unique(tokens["asym_id"])
            for asym_id in unique_chains:
                # Check if this is a protein chain
                chain_mask = tokens["asym_id"] == asym_id
                if np.any(chain_mask):
                    mol_type = tokens[chain_mask]["mol_type"][0]
                    if mol_type == const.chain_type_ids["PROTEIN"]:
                        _apply_mutations_to_tokens(tokens, chain_mask, chain_mutations_list, f"chain_{asym_id}")
        else:
            # Find chain by name
            chain_found = False
            for record_chain in tokenized.record.chains:
                if record_chain.chain_name == chain_id:
                    asym_id = record_chain.chain_id
                    chain_mask = tokens["asym_id"] == asym_id
                    if np.any(chain_mask):
                        _apply_mutations_to_tokens(tokens, chain_mask, chain_mutations_list, chain_id)
                        chain_found = True
                    break
            if not chain_found:
                raise ValueError(f"Chain {chain_id} not found in tokenized structure")
    
    # Create new tokenized object
    new_tokenized = Tokenized(
        tokens=tokens,
        bonds=tokenized.bonds,
        structure=tokenized.structure,
        msa=tokenized.msa,
        record=tokenized.record,
        residue_constraints=tokenized.residue_constraints,
        templates=tokenized.templates,
        template_tokens=tokenized.template_tokens,
        template_bonds=tokenized.template_bonds,
        extra_mols=tokenized.extra_mols,
    )
    
    return new_tokenized


def _apply_mutations_to_tokens(
    tokens: np.ndarray,
    chain_mask: np.ndarray,
    mutations: List[Tuple[int, str, str]],
    chain_name: str
) -> None:
    """Apply mutations to token array for a specific chain."""
    chain_tokens = tokens[chain_mask]
    chain_residue_indices = chain_tokens["residue_index"]
    
    for position, original_aa, new_aa in mutations:
        # Find token at this residue position
        token_mask = chain_residue_indices == position
        if not np.any(token_mask):
            raise ValueError(f"Position {position + 1} not found in chain {chain_name}")
            
        # Get new token ID
        new_token = const.prot_letter_to_token[new_aa]
        new_token_id = const.token_ids[new_token]
        
        # Apply mutation
        global_indices = np.where(chain_mask)[0][token_mask]
        for idx in global_indices:
            old_token_id = tokens[idx]["res_type"]
            tokens[idx]["res_type"] = new_token_id
            logging.info(f"Applied mutation to {chain_name}:{position + 1} (token {old_token_id} -> {new_token_id})")


def apply_masking_to_msa(
    msa_dict: Dict[str, MSA],
    mask_positions: List[int],
    mask_token: str = "X",
    mask_deletion_matrix: bool = True,
    chain_id: Optional[str] = None
) -> Dict[str, MSA]:
    """Apply masking to MSA features.
    
    Parameters
    ----------
    msa_dict : Dict[str, MSA]
        Dictionary of MSA objects by chain ID
    mask_positions : List[int]
        List of 0-based positions to mask
    mask_token : str
        Token to use for masking (default "X" for UNK)
    mask_deletion_matrix : bool
        Whether to zero deletion matrix at masked positions
    chain_id : Optional[str]
        Specific chain to mask, if None masks all chains
        
    Returns
    -------
    Dict[str, MSA]
        Dictionary with masked MSA objects
    """
    if not mask_positions:
        return msa_dict
        
    # Get mask token ID
    if mask_token not in const.prot_letter_to_token:
        raise ValueError(f"Invalid mask token: {mask_token}")
    mask_token_name = const.prot_letter_to_token[mask_token]
    mask_token_id = const.token_ids[mask_token_name]
    
    masked_msa_dict = {}
    
    for msa_chain_id, msa in msa_dict.items():
        if chain_id is not None and msa_chain_id != chain_id:
            # Skip this chain if specific chain requested
            masked_msa_dict[msa_chain_id] = msa
            continue
            
        # Create copies of MSA arrays
        new_residues = msa.residues.copy()
        new_deletions = msa.deletions.copy() if mask_deletion_matrix else msa.deletions
        new_sequences = msa.sequences.copy()
        
        # Get sequence length from first sequence
        if len(msa.sequences) == 0:
            masked_msa_dict[msa_chain_id] = msa
            continue
            
        first_seq = msa.sequences[0]
        seq_length = first_seq["res_end"] - first_seq["res_start"]
        
        # Validate mask positions
        valid_positions = [pos for pos in mask_positions if 0 <= pos < seq_length]
        if len(valid_positions) != len(mask_positions):
            invalid_positions = [pos + 1 for pos in mask_positions if pos < 0 or pos >= seq_length]
            logging.warning(
                f"Chain {msa_chain_id}: Skipping invalid positions {invalid_positions} "
                f"(sequence length: {seq_length})"
            )
        
        if not valid_positions:
            masked_msa_dict[msa_chain_id] = msa
            continue
            
        # Apply masking to non-target sequences (skip first sequence which is the query)
        for seq_idx in range(1, len(msa.sequences)):
            sequence = msa.sequences[seq_idx]
            res_start = sequence["res_start"]
            res_end = sequence["res_end"]
            
            # Mask specified positions
            for pos in valid_positions:
                if pos < (res_end - res_start):
                    global_res_idx = res_start + pos
                    if global_res_idx < len(new_residues):
                        new_residues[global_res_idx]["res_type"] = mask_token_id
        
        # Mask deletion matrix if requested
        if mask_deletion_matrix and len(valid_positions) > 0:
            # Filter out deletions at masked positions
            deletion_mask = np.ones(len(new_deletions), dtype=bool)
            for i, deletion in enumerate(msa.deletions):
                if deletion["res_idx"] in valid_positions:
                    deletion_mask[i] = False
            new_deletions = new_deletions[deletion_mask]
            
            # Update sequence deletion indices
            for seq_idx in range(len(new_sequences)):
                sequence = new_sequences[seq_idx]
                del_start = sequence["del_start"]
                del_end = sequence["del_end"]
                
                # Count how many deletions were removed before this sequence's deletions
                removed_before = np.sum(~deletion_mask[:del_start])
                removed_within = np.sum(~deletion_mask[del_start:del_end])
                
                new_sequences[seq_idx]["del_start"] = del_start - removed_before
                new_sequences[seq_idx]["del_end"] = del_end - removed_before - removed_within
        
        # Create new MSA object
        masked_msa = MSA(
            residues=new_residues,
            deletions=new_deletions,
            sequences=new_sequences,
        )
        masked_msa_dict[msa_chain_id] = masked_msa
        
        logging.info(f"Applied masking to chain {msa_chain_id} at positions {[p+1 for p in valid_positions]}")
    
    return masked_msa_dict


def validate_masking_config(masking_config: Dict) -> Dict:
    """Validate and normalize masking configuration.
    
    Parameters
    ----------
    masking_config : Dict
        Raw masking configuration
        
    Returns
    -------
    Dict
        Validated and normalized configuration
    """
    config = masking_config.copy()
    
    # Set defaults
    if "mask_token" not in config:
        config["mask_token"] = "X"
    if "mask_msa" not in config:
        config["mask_msa"] = True
    if "mask_deletion_matrix" not in config:
        config["mask_deletion_matrix"] = True
    if "mutations" not in config:
        config["mutations"] = []
        
    # Validate mask token
    if config["mask_token"] not in const.prot_letter_to_token:
        raise ValueError(f"Invalid mask token: {config['mask_token']}")
        
    return config 