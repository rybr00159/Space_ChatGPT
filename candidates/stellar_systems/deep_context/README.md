# Deeper independent context for the two stellar-system leads

Checked 2026-10-14 UTC. The [full context report](../../../reports/stellar_systems/deep_context.md)
documents properties, attempts to falsify interpretations, finite prior
recognition and positive query controls.

- **TIC 273879360 was previously publicly machine-flagged for eclipse-like
  features** in Kostov2025 Table2. It was unvetted and had no published period
  in that entry; this is not a first recognition of an entirely unflagged
  system. Gaia/infrared properties support a dwarf-like combined source but
  do not measure the stellar components. Its large faint-source RV error is
  not an epoch-RV variability detection.
- **TIC 260352037 is independently a spectroscopic red giant**, with Gaia
  GSP-Spec log g 2.51 and FLAME_spec radius 9.447 solar radii. Its measured RVs
  show no significant reported variation. A short-period detached system on
  the catalogued giant is implausible; the nearby faint stars remain source
  hypotheses, not spatially established hosts.
- **TIC 737466840 could supply enough light for the diluted eclipse**, but
  its colors are contaminated and compatible common-motion astrometry does
  not prove a bound wide pair. The central relative transverse velocity is
  above illustrative escape speeds, with errors too large to establish
  whether the pair is unbound.

[context_summary.json](context_summary.json) retains all numerical properties,
relative-motion covariance, conditional flux budgets, catalog controls and
42 verified raw-response source hashes. [source_properties.csv](source_properties.csv)
is a searchable concise table. Public raw queries and immutable manifests
are stored under `data/stellar_systems/deep_context/`.

The [independent ground-photometry analysis](ground_photometry/README.md)
checks earlier ASAS-SN and SkyMapper observations of source 273879360.

Reproduce retained-data derivation with:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python scripts/stellar_deep_context.py
```

Add `--download` to fetch missing public query snapshots. Raw failures are
retained and are never interpreted as absence of a catalog classification.
