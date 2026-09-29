# The method, parameter by parameter

This file covers the voxel descriptor and its tests. The steps before it (cells, folding, crystal benchmark, frames,
landmarks, the vector and arc30 reference methods) are described stage by stage in `upstream.md`.

## Structures and frame
- Input: one IMGT-numbered model per receptor, chains A (α) and B (β), from TCRBuilder2+ (ImmuneBuilder) with OpenMM
  refinement in upstream stage B. Heavy atoms only; altloc blank or `A`; hydrogens
  dropped. Receptors whose V-domain protein sequences are identical may be collapsed to one *molecule* (the clustering
  and label test then run on molecules).
- **Per-chain frame (F1):** each chain is Kabsch-fitted (Kabsch 1976) on its five framework landmarks — Cα of IMGT 23,
  41, 89, 104, 118 — onto the same chain of a fixed reference receptor, then placed on the TCR-intrinsic axes supplied in
  the frame file. Benchmark crystals and their models are placed with the same landmarks by their own parser.
- A whole-molecule frame (F0, ten landmarks) and framework-anchor variants exist for the frame diagnostic only.

## Atom sets (IMGT numbering; Lefranc et al. 2003)
- `loops`: CDR1 27–38, CDR2 56–65, HV4 81–86, CDR3 105–117, both chains, insertion codes included.
- `cdr3`: 105–117 only. Cys104 and Phe/Trp118 count as framework.

## Grid
- Voxel edge h = 1.0 Å. Per-chain boxes are frozen before the run: for loop arms, a margin (default 7 Å = 3σ + 1 Å)
  around the repertoire envelope; for CDR3 arms, the full min/max envelope of every CDR3 atom in the repertoire **and**
  the benchmark, plus the same margin, so no atom is truncated. Check: every atom ≥ 2σ from every box face.
- Each atom is an isotropic Gaussian of width σ (2.0 Å default; 1.5 Å in the sensitivity arms), point-sampled at voxel
  centres, truncated per axis at 3σ, with each 1-D factor renormalised to sum 1, so an atom deposits exactly its channel
  mass unless the box clips it. Sampling adequacy: keep σ/h ≥ 1.5; a smaller σ needs a smaller voxel.
- Channels and mass per atom (the typing table in `vxgrid.py`), **unweighted in the distance**:

| channel | atoms | mass |
|---|---|---|
| occupancy | every heavy atom | 1 |
| hydrophobic | carbons not bonded to N or O | 1 |
| aromatic | ring atoms of Phe, Tyr, Trp, His | 1 |
| H-bond donor | backbone N except Pro; Arg NE/NH1/NH2, Asn ND2, Gln NE2, Lys NZ, Trp NE1, His ND1/NE2, Ser OG, Thr OG1, Tyr OH | 1 |
| H-bond acceptor | every O; His ND1/NE2 | 1 |
| positive | Lys NZ; Arg CZ/NE/NH1/NH2; His ND1/NE2 | 1 / 0.25 each / 0.05 each |
| negative | Asp OD1/OD2, Glu OE1/OE2; C-terminal OXT | 0.5 each / 1 |

  A side-chain channel is stored but never used in the distance. Because nothing is normalised per channel, a channel's
  influence follows how many atoms carry it and how their mass is distributed: charge spread over four atoms contributes
  less than the same charge on one atom, and occupancy usually dominates. The voxelised pharmacophore scheme follows
  Jiménez et al. 2017 (DeepSite).
- Grids are stored float16 (row = receptor) and deleted after their last reader (`keep_grids` keeps them; `pipeline.md`).

## Distance
- D² = Σ over the two chain boxes, Σ over channels, Σ over voxels (a − b)² / h³ — computed as (|a|² + |b|² − 2a·b)/h³ with
  the Gram matrix accumulated over feature blocks in float64.
- It splits exactly into 14 non-negative parts (2 chains × 7 channels); the report shows each member's mean over the rest
  of its cluster, and a gate checks that the parts sum to D².
- Build checks, all stop-on-failure: typing complete; mass conservation; grid distance vs the exact grid-free (analytic
  Gaussian) distance on random pairs (< 2 %); fast grid equals a slow reference implementation; permuting atom identities
  leaves occupancy unchanged and changes every chemistry channel; blocked Gram equals `pdist`; control pairs exact;
  triangle inequality.

## Threshold and clustering
- Draw once, with a fixed seed, and reuse across arms: a set of near-identical control pairs (one CDR3 residue apart), a
  large random background sample (default 60,000 pairs), and a within-CDR3-length-class sample.
- Threshold = 1st percentile of the background distances, rounded to 4 dp. Complete linkage on the molecule distance
  matrix; clusters are numbered by size; singletons dropped. Stability: adjusted Rand index (Hubert & Arabie 1985) at
  threshold × 0.98 and × 1.02.
- Reported, not gates (unless you pre-register them as gates): the share of control pairs within the threshold, and the
  crystal-model error over the threshold.

## Validity panel (crystal benchmark)
- Receptors with both a crystal structure and a model of the same sequence; comparisons use only atoms present in both,
  with per-pair masks for positions missing in either structure.
- C = crystal–crystal distances, M = model–model distances, S = a sequence baseline (BLOSUM62 per loop; Henikoff &
  Henikoff 1992).
- **E1** = partial Spearman of C and M given S (structure recovered from models beyond sequence). **E2** = the same given
  a reference descriptor's M and S (structure the reference misses). **E4** = 1 − 2·mean d²(crystal, own model) / mean
  d²(background pairs) (reliability). Plus the median error over the pilot threshold.
- Intervals: bootstrap over receptors (default 2,000 resamples, fixed seed). For paired comparisons (for example two
  blurs) use the **same** resamples for both and report the difference distribution — it is tighter than comparing two
  intervals by eye. Mantel's test (Mantel 1967) for plain matrix correlations.

## Label test
- Runs once per arm, after clusters are fixed. On molecules with a label:
  - **Purity** over clusters ≥ 3, against a null ladder of permutations (default 200) with a fixed seed: unstratified,
    within donor, within donor + V-gene pair, and within donor + V pair + CDR3 length class. Report excess and z per rung
    — the V-matched rung is the one that says whether shape adds anything beyond germline.
  - **Per cluster** (≥ 2 members): one-sided binomial of the most common label against its repertoire base rate;
    Benjamini–Hochberg (1995) across clusters within the arm; survivors at q ≤ 0.15 (set your own threshold in advance).
  - **Overlap** with a reference method's clusters by a match rule (default: ≥ 75 % of members in one cluster and
    significant), in both directions.
- Validate the implementation first by reproducing a reference method's published cluster tests exactly.

## Sequence cross-check (tcrdist3)
- Paired αβ TCRdist (Dash et al. 2017; tcrdist3, Mayer-Blackwell et al. 2021) with package defaults: CDR1, CDR2, CDR2.5
  weight 1 with a fixed gap position; CDR3 weight 3, trimmed 3 residues N-terminal and 2 C-terminal, best gap position;
  gap penalty 4; α + β summed. Genes are mapped to the package's own database (as `*01` alleles); genes it lacks are
  counted and excluded, never guessed. In tcrdist3 0.3 the IMGT gap symbol in germline CDR strings scores 0 against any
  residue. Above ~10,000 clones the dense matrix must be requested explicitly.
- An independent partition is built with the **same** threshold rule as the voxel clusters, so the only difference is the
  distance. Each structural cluster is then *also* (≥ 75 % of members in one sequence cluster), *partly* (50–75 %) or
  *not* a sequence cluster, plus descriptive pairwise statistics.
- Verify tcrdist3 against an independent implementation of the published formula on random pairs, and check that the
  near-identical control pairs fall far below the threshold.

## Reproducibility conventions
- Fix every random seed in the config or the script; draw pair samples once and reuse them across arms.
- Never change a finished stage's logic. Shared helpers may gain optional arguments whose defaults keep earlier calls
  identical (the blur argument works this way).
- Store distance matrices, not grids; every grid is regenerable from its script.

## References
Kabsch W. *Acta Crystallogr A* 1976;32:922–3 · Lefranc M-P, et al. *Dev Comp Immunol* 2003;27:55–77 (IMGT) ·
Jiménez J, et al. DeepSite. *Bioinformatics* 2017;33:3036–42 (PMID 28575181) · Ragoza M, et al. *J Chem Inf Model*
2017;57:942–57 (PMID 28368587) · Grant JA, Pickup BT. *J Phys Chem* 1995;99:3503–10 (not PubMed-indexed) ·
Henikoff S, Henikoff JG. *PNAS* 1992;89:10915–9 · Mantel N. *Cancer Res* 1967;27:209–20 · Hubert L, Arabie P.
*J Classif* 1985;2:193–218 · Benjamini Y, Hochberg Y. *J R Stat Soc B* 1995;57:289–300 · Koo TK, Li MY. *J Chiropr Med*
2016;15:155–63 · Dash P, et al. *Nature* 2017;547:89–93 (PMID 28636592) · Mayer-Blackwell K, et al. *eLife*
2021;10:e68605 (PMID 34845983) · Dunbar J, et al. *PLoS Comput Biol* 2014;10:e1003852 (TRangle) · Quast NP, et al.
*Commun Biol* 2025;8:362 (TCRBuilder2+) · Abanades B, et al. *Commun Biol* 2023;6:575 (ImmuneBuilder) · Rego N, Koes D. *Bioinformatics* 2015;31:1322–4 (PMID 25505090, 3Dmol.js).
