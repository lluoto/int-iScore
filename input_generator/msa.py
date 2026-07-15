"""
MSA computation module.

Connects to ColabFold MMseqs2 server for multiple sequence alignment.
Supports AF3, Boltz, and Chai MSA file formats.
"""
import os
import time
import logging

logger = logging.getLogger(__name__)

try:
    from mmseqs2 import run_mmseqs2
except ImportError:
    run_mmseqs2 = None


def parse_a3m_sequences(a3m_text):
    """Parse a3m text into (headers, sequences) lists."""
    lines = a3m_text.strip().splitlines()
    headers, sequences = [], []
    current_header = None
    current_parts = []
    for line in lines:
        if line.startswith('>'):
            if current_header is not None:
                sequences.append(''.join(current_parts))
                current_parts = []
            headers.append(line)
            current_header = line
        else:
            current_parts.append(line.strip())
    if current_header is not None:
        sequences.append(''.join(current_parts))
    return headers, sequences


def compute_msa(sequence, target_id, msa_dir, host_url='https://api.colabfold.com',
                use_env=True, use_filter=True, force=False):
    """
    Compute MSA for a protein sequence via ColabFold MMseqs2 server.

    Args:
        sequence: Amino acid sequence string
        target_id: Unique identifier for the protein
        msa_dir: Output directory for MSA files
        host_url: ColabFold server URL
        use_env: Include environmental MSAs
        use_filter: Filter MSAs
        force: Recompute even if cached

    Returns:
        (a3m_path, a3m_text) or (None, None) on failure
    """
    if run_mmseqs2 is None:
        raise ImportError('mmseqs2 module not found. Install from ColabFold.')

    os.makedirs(msa_dir, exist_ok=True)
    a3m_path = os.path.join(msa_dir, '{}.a3m'.format(target_id))

    if os.path.exists(a3m_path) and not force:
        with open(a3m_path) as f:
            return a3m_path, f.read()

    tmp_prefix = os.path.join(msa_dir, '{}_tmp'.format(target_id))
    try:
        a3m_lines, _ = run_mmseqs2(
            sequence, prefix=tmp_prefix,
            use_env=use_env, use_filter=use_filter,
            use_pairing=False, host_url=host_url,
        )
        a3m_text = a3m_lines[0]
        with open(a3m_path, 'w') as f:
            f.write(a3m_text)
        # Cleanup temp
        for p in ['{}_env'.format(tmp_prefix), '{}_all'.format(tmp_prefix)]:
            if os.path.exists(p):
                import shutil
                shutil.rmtree(p, ignore_errors=True)
        return a3m_path, a3m_text
    except Exception as e:
        logger.error('MSA failed for {}: {}'.format(target_id, e))
        return None, None


def batch_compute_msa(sequences, msa_dir, host_url='https://api.colabfold.com',
                      workers=1, force=False):
    """
    Compute MSAs for multiple sequences.

    Args:
        sequences: dict of {target_id: sequence}
        msa_dir: Output directory
        host_url: ColabFold server URL
        workers: Parallel workers (1=sequential)
        force: Recompute cached

    Returns:
        dict of {target_id: (a3m_path, a3m_text)}
    """
    results = {}
    for idx, (tid, seq) in enumerate(sorted(sequences.items())):
        logger.info('MSA [{}/{}] {} ({} aa)'.format(idx + 1, len(sequences), tid, len(seq)))
        path, text = compute_msa(seq, tid, msa_dir, host_url, force=force)
        results[tid] = (path, text)
        if path:
            time.sleep(1)  # rate limiting
    return results
