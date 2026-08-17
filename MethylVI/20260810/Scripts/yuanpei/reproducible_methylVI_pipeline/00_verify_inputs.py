#!/usr/bin/env python3
import json
from pathlib import Path
import re

import anndata as ad
import pandas as pd

DATA = Path("/home/lijia/jiangyuanpei/methscan/xunyin/20260409_mix_0513/allcools_5kbin")
OUT = DATA / "methylVI/reproducible_5kbin_pipeline/input_audit.json"


def canonical(x):
    return re.sub(r"^([^-_]+)-", r"\1_", str(x))


def main():
    h5ad = DATA / "mcg_5kb.clustered.h5ad"
    anno_path = DATA / "cell_type_manual_with_donor_disease.csv"
    allc_dir = DATA / "input_allc"
    a = ad.read_h5ad(h5ad, backed="r")
    anno = pd.read_csv(anno_path, dtype=str)
    required = ["cell_id", "cell_type", "donor", "disease"]
    if list(anno.columns[:4]) != required:
        raise ValueError(f"Annotation columns must start with {required}")
    anno["match_id"] = anno.cell_id.map(canonical)
    if anno.match_id.duplicated().any():
        raise ValueError("Duplicated annotation IDs after '-'/'_' normalization")
    amap = anno.set_index("match_id")
    cells = pd.Index(a.obs_names.astype(str))
    normalized = cells.map(canonical)
    matched = normalized.isin(amap.index)
    prefixes = normalized.str.extract(r"^([^_]+)_", expand=False)
    valid_donors = set(anno.donor.dropna())
    donor_assignable = matched | prefixes.isin(valid_donors)
    files = list(allc_dir.glob("*.allc.tsv.gz"))
    file_ids = {canonical(p.name.removesuffix(".allc.tsv.gz")) for p in files}
    allc_matched = normalized.isin(file_ids)
    audit = {
        "h5ad": str(h5ad), "cells": int(a.n_obs), "retained_5kb_bins": int(a.n_vars),
        "h5ad_layers": list(a.layers.keys()), "latent_source": "AllCools filtered 5kb bins",
        "annotation": str(anno_path), "annotation_rows": len(anno),
        "fully_annotated_selected_cells": int(matched.sum()),
        "donor_prefix_inferred_cells": int((~matched & donor_assignable).sum()),
        "donor_assignable_cells": int(donor_assignable.sum()),
        "allc_files": len(files), "selected_cells_with_allc": int(allc_matched.sum()),
        "id_normalization": "replace only the first donor/barcode '-' delimiter with '_'",
    }
    if not donor_assignable.all() or not allc_matched.all():
        raise ValueError(f"Input audit failed: {audit}")
    OUT.write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
