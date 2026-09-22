#!/usr/bin/env python3
"""Read Kestrel NSW=0 complexes and relaxed references; audit new SPE energies."""
import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal
import functools
import gzip
import hashlib
import json
import math
from pathlib import Path
import re

from ase.io import read
from bs4 import BeautifulSoup
import numpy as np

from add_missing_dft_images import incar_values, normalize_setting
from extract_perlmutter_energies import parse_output
from molecule_names import canonical

ROOT = Path(__file__).resolve().parents[1]
FUNC_DIRS = {"PBE": "pbe", "PBE_D3": "pbe_d3", "r2scan": "r2scan", "beef_vdw": "beef_vdw"}
MAX_ADS_EV = 5.0


def infer_functional(settings):
    if normalize_setting("METAGGA", settings.get("METAGGA", "")) == "R2SCAN":
        return "r2scan"
    if settings.get("GGA", "PE").upper() == "BF" and normalize_setting("LUSE_VDW", settings.get("LUSE_VDW", "")) == "T":
        return "beef_vdw"
    if settings.get("GGA", "PE").upper() == "PE" and not settings.get("METAGGA"):
        if normalize_setting("LUSE_VDW", settings.get("LUSE_VDW", "F")) == "T":
            return "unknown"
        if settings.get("IVDW", "0") in ("11", "12"):
            return "pbe_d3"
        if settings.get("IVDW", "0") == "0":
            return "pbe"
    return "unknown"


@functools.lru_cache(None)
def component(directory):
    result = {"directory": str(directory), "incar": str(directory / "INCAR"),
              "outcar": str(directory / "OUTCAR"), "status": "missing", "energy": None,
              "nsw": None, "composition": {}, "potentials": {}, "functional": "unknown"}
    try:
        settings = incar_values((directory / "INCAR").read_text())
        result.update(settings=settings, functional=infer_functional(settings), nsw=int(settings["NSW"]))
        path = directory / "OUTCAR"
        if not path.exists():
            path = directory / "OUTCAR.gz"
        result["outcar"] = str(path)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rb") as source:
            head_bytes = source.read(100000)
            if path.suffix == ".gz":
                tail_bytes = head_bytes
                for chunk in iter(lambda: source.read(1000000), b""):
                    tail_bytes = (tail_bytes + chunk)[-1000000:]
            else:
                source.seek(max(0, path.stat().st_size - 1000000))
                tail_bytes = source.read()
        head, tail = head_bytes.decode(errors="replace"), tail_bytes.decode(errors="replace")
        parsed = parse_output(head, tail)
        result.update(energy=parsed["energy"], status=parsed["convergence"], composition=parsed["composition"],
                      potentials=parsed["potentials"], outcar_size=path.stat().st_size,
                      outcar_mtime_ns=path.stat().st_mtime_ns,
                      outcar_head_tail_sha256=hashlib.sha256(head_bytes + b'\0' + tail_bytes).hexdigest())
        if parsed.get("failure_reason"):
            result["note"] = parsed["failure_reason"]
        for key in ("NSW", "IBRION", "EDIFF", "EDIFFG", "GGA", "IVDW", "METAGGA", "LUSE_VDW"):
            if key not in settings:
                continue
            values = re.findall(r"\b" + key + r"\s*=\s*([^\s;]+)", head)
            if key in ("NSW", "IBRION", "EDIFF", "EDIFFG", "IVDW"):
                matches = bool(values) and math.isclose(float(values[-1]), float(settings[key]), rel_tol=1e-6, abs_tol=1e-12)
            else:
                matches = bool(values) and normalize_setting(key, values[-1]) == normalize_setting(key, settings[key])
            if not matches:
                raise ValueError(f"INCAR/OUTCAR {key} mismatch")
        if not parsed["metadata_ok"]:
            raise ValueError("Missing OUTCAR composition/potential metadata")
        if parsed["energy"] is None or not math.isfinite(parsed["energy"]):
            raise ValueError("Missing/nonfinite final TOTEN")
        atoms = read(directory / "POSCAR", format="vasp")
        if dict(Counter(atoms.get_chemical_symbols())) != parsed["composition"]:
            raise ValueError("POSCAR/OUTCAR composition mismatch")
        if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell).all():
            raise ValueError("Nonfinite POSCAR geometry")
        result["cell"] = atoms.cell.tolist()
        if result["functional"] == "unknown":
            raise ValueError("Unrecognized functional")
    except (OSError, ValueError, KeyError, IndexError) as error:
        result.update(status="invalid", note=str(error))
    from structure_identity import attach_geometry
    return attach_geometry(result,directory)


def calculation_dirs(system):
    for base in (system, system / "singlepoint"):
        if (base / "INCAR").is_file():
            yield base, None
        for name, functional in FUNC_DIRS.items():
            directory = base / name
            if (directory / "INCAR").is_file():
                yield directory, functional


def is_relaxed(result, functional):
    return (result["status"] == "converged" and result["functional"] == functional
            and result["nsw"] is not None and result["nsw"] > 0
            and result.get("settings", {}).get("IBRION") in ("1", "2", "3"))


def assess(complex_result, slab, molecule, functional):
    if complex_result["nsw"] != 0:
        return "not_singlepoint", "Complex NSW is not zero", None
    if complex_result["functional"] != functional:
        return "functional_mismatch", "Complex functional disagrees with directory", None
    if complex_result["status"] != "converged":
        return "complex_unconverged", complex_result.get("note", complex_result["status"]), None
    if slab is None or molecule is None:
        return "reference_unavailable", "No converged relaxed references with matching functional/composition", None
    if not is_relaxed(slab, functional) or not is_relaxed(molecule, functional):
        return "reference_unavailable", "References must be converged relaxations of the same functional", None
    expected = Counter(slab["composition"]) + Counter(molecule["composition"])
    if dict(expected) != complex_result["composition"]:
        return "composition_mismatch", "Complex does not equal slab plus molecule composition", None
    energy = float(Decimal(str(complex_result["energy"])) - Decimal(str(slab["energy"])) - Decimal(str(molecule["energy"])))
    from potential_matching import potential_match
    if not all(potential_match(complex_result,r) for r in [slab,molecule]):
        return 'potential_mismatch','Complex and reference PAW potential identities differ or are missing; energy retained only as diagnostic',energy
    if abs(energy) > MAX_ADS_EV:
        return "energy_review", "|E_ads| > 5 eV; recorded for review, not a convergence failure", energy
    return "ok", "NSW=0 complex; converged relaxed references; matching composition", energy


def write_csv(path, rows, fields):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output-root", type=Path, default=ROOT)
    args = parser.parse_args()
    project, output = args.project_root.resolve(), args.output_root.resolve()
    page = BeautifulSoup((output / "dft_comparison.html").read_text(), "html.parser")
    targets = {(card["data-surf"], canonical(card["data-mol"])) for card in page.select(".g")}
    surfaces, molecules = {pair[0] for pair in targets}, {pair[1] for pair in targets}
    references = defaultdict(list)
    components = {}
    for role, roots in (("molecule", [project / "vasp_mol"]),
                        ("slab", [project / "vasp_slab_kestrel", project / "vasp_slab"])):
        for root in roots:
            if not root.is_dir():
                continue
            for system in sorted(root.iterdir()):
                if not system.is_dir():
                    continue
                label = canonical(system.name) if role == "molecule" else re.sub(r"_n\d+$", "", system.name)
                if label not in (molecules if role == "molecule" else surfaces):
                    continue
                for directory, expected_func in calculation_dirs(system):
                    result = component(directory)
                    components[str(directory)] = result
                    if expected_func is not None and result["functional"] != expected_func:
                        continue
                    system_label = result.get("settings", {}).get("SYSTEM", "")
                    if not (canonical(system_label) == label if role == "molecule" else system_label in (system.name, label)):
                        continue
                    if is_relaxed(result, result["functional"]):
                        references[(role, label, result["functional"])].append(result)
    print(f"Audited {len(components)} reference calculations; {sum(map(len, references.values()))} converged relaxed references.", flush=True)
    best = project / "poscar" / "best"
    if not best.is_dir():
        raise SystemExit("No poscar/best directory; existing exports were preserved. Use component_store.py for other layouts.")
    systems = []
    for entry in sorted(best.iterdir()):
        if entry.is_dir() and re.fullmatch(r"C\d+", entry.name):
            systems.extend(sorted(path for path in entry.iterdir() if path.is_dir()))
        elif entry.is_dir():
            systems.append(entry)
    rows = []
    last_reported = 0
    for system in systems:
        pair = system.name.split("_", 1)
        if len(pair) != 2:
            continue
        surface, raw_molecule = pair
        molecule = canonical(raw_molecule)
        if (surface, molecule) not in targets:
            continue
        metal = re.match(r"[A-Z][a-z]?", surface).group()
        for directory, expected_func in calculation_dirs(system):
            comp = component(directory)
            components[str(directory)] = comp
            functional = expected_func or comp["functional"]
            if functional not in FUNC_DIRS.values():
                continue
            slab_candidates = [r for r in references[("slab", surface, functional)]
                               if r["composition"] == {metal: comp["composition"].get(metal, -1)}]
            adsorbate = dict(comp["composition"])
            adsorbate.pop(metal, None)
            mol_candidates = [r for r in references[("molecule", molecule, functional)] if r["composition"] == adsorbate]
            # Prefer the same in-plane cell, then the explicitly requested reference root.
            slab_candidates.sort(key=lambda r: (not np.allclose(np.array(r.get("cell", np.zeros((3,3))))[:2], np.array(comp.get("cell", np.ones((3,3))))[:2], atol=1e-5, rtol=0),
                                                "vasp_slab_kestrel" not in r["directory"], r["directory"]))
            mol_candidates.sort(key=lambda r: (Path(r["directory"]).parent.name != raw_molecule, r["directory"]))
            slab = slab_candidates[0] if slab_candidates else None
            mol = mol_candidates[0] if mol_candidates else None
            status, note, energy = assess(comp, slab, mol, functional)
            if comp.get("settings", {}).get("SYSTEM") != system.name:
                status, note = "system_mismatch", "INCAR SYSTEM disagrees with exact system directory"
            row = {"surface": surface, "molecule": molecule, "functional": functional,
                   "status": status, "note": note, "E_ads_SPE": energy,
                   "complex_directory": str(directory), "complex_NSW": comp["nsw"],
                   "slab_directory": slab["directory"] if slab else "",
                   "molecule_directory": mol["directory"] if mol else "",
                   "E_complex": comp["energy"], "E_slab_relaxed": slab["energy"] if slab else None,
                   "E_molecule_relaxed": mol["energy"] if mol else None,
                   "slab_NSW": slab["nsw"] if slab else None, "molecule_NSW": mol["nsw"] if mol else None,
                   "matching_slab_references": len(slab_candidates), "matching_molecule_references": len(mol_candidates),
                   "selected": "false"}
            rows.append(row)
        if len(rows) >= last_reported + 200:
            print(f"Audited {len(rows)} complex candidates...", flush=True)
            last_reported = len(rows)
    selected = {}
    for row in rows:
        if row["status"] != "ok":
            continue
        key = (row["surface"], row["molecule"], row["functional"])
        def priority(candidate):
            path = Path(candidate["complex_directory"])
            return ("singlepoint" not in path.parts, len(path.parts), str(path))
        if key not in selected or priority(row) < priority(selected[key]):
            selected[key] = row
    for row in selected.values():
        row["selected"] = "true"
    if not rows:
        raise SystemExit("No matching candidates; existing exports were preserved.")
    fields = list(rows[0])
    write_csv(output / "dft_kestrel_singlepoint_audit.csv", rows, fields)
    write_csv(output / "dft_kestrel_singlepoint.csv", [selected[key] for key in sorted(selected)], fields)
    (output / "dft_kestrel_singlepoint_sources.json").write_text(json.dumps({
        "formula": "E_ads_SPE = E_complex_NSW0 - E_slab_relaxed - E_molecule_relaxed",
        "policy": "Preserve existing energies; only fill missing SPE cells. Require completion, final electronic convergence, NSW=0 complexes, NSW>0 converged ionic references, functional and composition matches. Same |E_ads| <= 5 eV publication screen. Matching PAW TITEL identities are required; mismatched subtractions are diagnostic only. No energy scaling, offsets, gas-energy overrides, or sign-based rejection.",
        "complex_selection": "Prefer singlepoint directories, then shallowest path, then lexical path; do not select by energy.",
        "reference_selection": "Same surface/molecule, functional and atom counts. Slab: prefer matching in-plane cell, then vasp_slab_kestrel, then lexical path. Molecule: prefer exact original name, then lexical path. Only completed relaxed references qualify.",
        "components": components}, indent=2) + "\n")
    print(f"Audited {len(rows)} complex candidates: {dict(Counter(row['status'] for row in rows))}; {len(selected)} screened SPE results.")


if __name__ == "__main__":
    main()
