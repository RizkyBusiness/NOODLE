"""Two runtime fixes to ImmuneBuilder 1.x that any TCR repertoire will hit.

1. TRAV/DV genes (TRAV14/DV4, TRAV29/DV5, TRAV38-2/DV8, ...) are shared between
   the TRA and TRD loci, so ANARCI assigns those alpha V domains to the delta
   chain type. ImmuneBuilder's sequence check restricts numbering to the chain
   type it was handed and therefore rejects genuine alpha chains using these
   genes - about 15% of clonotypes in the repertoire this pipeline came from.
   We widen the allowed set for "A" to include "D" and vice versa, mirroring the
   L/K handling already present for antibody light chains.

2. ImmuneBuilder.refine.strained_sidechain_bonds_fixer builds the OpenMM
   platform properties as a set literal, {'Threads', str(n_threads)}, instead of
   a dict. Any model that reaches the strained-bond fixer with n_threads > 0
   crashes in Context construction. We rebuild the function with the dict.

Both patches are checked before they are applied and are no-ops on versions
that no longer need them, so this stays safe as ImmuneBuilder moves on.
"""
import inspect


def apply():
    """Patch in place; returns the list of fixes actually applied."""
    fixed = []

    import ImmuneBuilder.sequence_checks as sc
    if getattr(sc, "_patched_here", False):
        return ["already_applied"]
    src = inspect.getsource(sc.number_single_sequence)
    old = '    if chain == "L":\n        allow.append("K")\n'
    new = old + ('    if chain == "A":\n        allow.append("D")\n'
                 '    if chain == "D":\n        allow.append("A")\n')
    if old in src and 'allow.append("D")' not in src:
        exec(compile(src.replace(old, new), "<ib_patch_chain>", "exec"), sc.__dict__)
        fixed.append("shared_locus_chain_types")

    import ImmuneBuilder.refine as rf
    src = inspect.getsource(rf.strained_sidechain_bonds_fixer)
    if "{'Threads', str(n_threads)}" in src:
        exec(compile(src.replace("{'Threads', str(n_threads)}",
                                 "{'Threads': str(n_threads)}"),
                     "<ib_patch_threads>", "exec"), rf.__dict__)
        fixed.append("openmm_threads_dict")

    sc._patched_here = True
    return fixed
