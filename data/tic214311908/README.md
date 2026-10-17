# Preserved observations and catalogue checks for TIC 214311908

Read the scientific assessment in `../../candidates/tic214311908/README.md`.

`observations.json` identifies four original QLP products. Catalogue, resolved-image, official-method and TARS recognition folders retain exact responses, requests and immutable hashes. The parent TARS catalogue is about136MB and remains outside ordinary Git, with its public URL, publisher MD5, SHA256 and exact target-scan result recorded. No downloaded dataset was discarded.

`checkpoint/` embeds **all13scientific TESS FITS** used in this review: four QLP light curves, four native TESSCut cubes (S66/93/103/104), four TARS curves and one ELEANOR-LITE curve. Bounded archives separate the native cubes and the nine light curves. A further archive preserves compact cache queries/metadata. The original redundant download ZIPs and parent catalogue remain cached and inventoried rather than duplicated in Git.

```sh
python scripts/tic214311908_checkpoint.py verify
python scripts/tic214311908_checkpoint.py restore --destination /tmp/tic214311908_new_offline
```

For a lost cache, restore into `/workspace/tess-research-data/tic214311908`. Then `verify --science-only` verifies all scientific FITS, sources and derived artifacts without requiring redundant ZIPs or the unneeded parent catalogue. Ordinary `verify` additionally checks every retained cached file. Existing mismatching files are never overwritten and expected input hashes are never replaced to accept changed bytes.

A fresh task actually started without the prior local caches, so this work does not assume cloud filesystem durability. Important scientific FITS restore byte for byte without network. Dynamic TESSCut re-creation need not reproduce original container/header bytes, which is why the actual cubes are archived directly. Source FFI identifiers, explicit timing standards, quality flags and original processing metadata remain preserved.
