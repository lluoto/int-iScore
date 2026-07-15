"""
Model input generators.

Generate input JSONs for AlphaFold3, Boltz-1, and Chai-1.
"""
import json
import random
import os


def make_af3_json(name, seq1, seq2=None, seeds=20, seed_base=42,
                  msa1=None, msa2=None, dialect='alphafoldserver'):
    """
    Create AlphaFold3 job_request JSON.

    Args:
        name: Job name
        seq1: First chain sequence
        seq2: Second chain sequence (None for monomer, same seq for homodimer)
        seeds: Number of random seeds
        seed_base: Base seed for RNG
        msa1: Unpaired MSA text for chain 1 (optional)
        msa2: Unpaired MSA text for chain 2 (optional)
        dialect: AF3 dialect

    Returns:
        List with single job_request dict
    """
    sequences = []

    # Chain 1
    chain1 = {'proteinChain': {'sequence': seq1, 'count': 1}}
    if msa1:
        chain1['proteinChain']['unpairedMsa'] = msa1

    if seq2 is None:
        sequences.append(chain1)
    elif seq2 == seq1:
        # Homodimer: single entry with count=2
        chain1['proteinChain']['count'] = 2
        if msa1:
            chain1['proteinChain']['unpairedMsa'] = msa1
        sequences.append(chain1)
    else:
        sequences.append(chain1)
        chain2 = {'proteinChain': {'sequence': seq2, 'count': 1}}
        if msa2:
            chain2['proteinChain']['unpairedMsa'] = msa2
        sequences.append(chain2)

    rng = random.Random(seed_base)
    seed_list = [str(rng.randint(1, 999999999)) for _ in range(seeds)]

    return [{
        'name': name,
        'modelSeeds': seed_list,
        'sequences': sequences,
        'dialect': dialect,
        'version': 1,
    }]


def make_boltz_input(name, seq1, seq2=None, msa_path1=None, msa_path2=None):
    """
    Create Boltz-1 input dict.
    """
    chains = []
    if seq2 is None:
        chains.append({'protein': {'id': ['A'], 'sequence': seq1,
                                    'msa': msa_path1 or ''}})
    elif seq2 == seq1:
        chains.append({'protein': {'id': ['A', 'B'], 'sequence': seq1,
                                    'msa': msa_path1 or ''}})
    else:
        chains.append({'protein': {'id': ['A'], 'sequence': seq1,
                                    'msa': msa_path1 or ''}})
        chains.append({'protein': {'id': ['B'], 'sequence': seq2,
                                    'msa': msa_path2 or ''}})

    return {
        'version': 2, 'name': name, 'sequences': chains,
        'seed': 42, 'recycling_steps': 3, 'sampling_steps': 200,
    }


def make_chai_input(name, seq1, seq2=None, msa_path1=None, msa_path2=None):
    """
    Create Chai-1 input dict.
    """
    sequences = []
    if seq2 is None:
        sequences.append({'protein': {'id': ['A'], 'sequence': seq1,
                                       'msa': msa_path1 or ''}})
    elif seq2 == seq1:
        sequences.append({'protein': {'id': ['A', 'B'], 'sequence': seq1,
                                       'msa': msa_path1 or ''}})
    else:
        sequences.append({'protein': {'id': ['A'], 'sequence': seq1,
                                       'msa': msa_path1 or ''}})
        sequences.append({'protein': {'id': ['B'], 'sequence': seq2,
                                       'msa': msa_path2 or ''}})

    return {'name': name, 'sequences': sequences, 'modelSeeds': [42], 'version': 2}


def generate_all_formats(name, seq1, seq2=None, output_dir='.',
                         msa_text1=None, msa_text2=None, msa_path1=None, msa_path2=None,
                         seeds=20, seed_base=42):
    """
    Generate all model input formats (AF3, Boltz, Chai) for one pair.

    Args:
        name: Base name for files
        seq1/seq2: Sequences
        output_dir: Output directory
        msa_text1/msa_text2: MSA text for AF3 JSON embedding
        msa_path1/msa_path2: MSA file paths for Boltz/Chai
        seeds: Number of AF3 seeds
        seed_base: Base seed

    Returns:
        dict of {format: filepath}
    """
    os.makedirs(output_dir, exist_ok=True)

    # AF3
    af3 = make_af3_json(name, seq1, seq2, seeds, seed_base, msa_text1, msa_text2)
    af3_path = os.path.join(output_dir, '{}_af3.json'.format(name))
    with open(af3_path, 'w') as f:
        json.dump(af3, f, indent=2)

    # Boltz
    boltz = make_boltz_input(name, seq1, seq2, msa_path1, msa_path2)
    boltz_path = os.path.join(output_dir, '{}_boltz.json'.format(name))
    with open(boltz_path, 'w') as f:
        json.dump(boltz, f, indent=2)

    # Chai
    chai = make_chai_input(name, seq1, seq2, msa_path1, msa_path2)
    chai_path = os.path.join(output_dir, '{}_chai.json'.format(name))
    with open(chai_path, 'w') as f:
        json.dump(chai, f, indent=2)

    return {'af3': af3_path, 'boltz': boltz_path, 'chai': chai_path}
