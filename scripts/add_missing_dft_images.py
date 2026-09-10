#!/usr/bin/env python3
"""Audit missing structure images using read-only access to Kestrel relaxations.

Requires ASE, NumPy, Matplotlib and BeautifulSoup4. Run without --apply to audit.
Only --apply writes website images, provenance and targeted HTML replacements.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import html
import json
import os
from pathlib import Path
import re

import ase
from ase.io import read, write
from ase.geometry import find_mic
from bs4 import BeautifulSoup
import numpy as np

from extract_perlmutter_energies import read_output

ROOT = Path(__file__).resolve().parents[1]
FUNCTIONALS = {"PBE": "PBE", "PBE_D3": "PBE+D3", "r2scan": "r2SCAN", "beef_vdw": "BEEF-vdW"}
FORMULAS = {"H2": {"H": 2}, "O2": {"O": 2}, "formate": {"C": 1, "H": 1, "O": 2},
            "formic_acid": {"C": 1, "H": 2, "O": 2}, "CO": {"C": 1, "O": 1},
            "methanol": {"C": 1, "H": 4, "O": 1}}
ROTATION = "-70x,20y,10z"


def incar_values(text):
    values = {}
    for line in text.splitlines():
        for setting in re.split(r"[!#]", line, maxsplit=1)[0].split(";"):
            if "=" in setting:
                key, value = setting.split("=", 1)
                values[key.strip().upper()] = value.strip()
    return values


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def final_positions_forces(tail):
    blocks = re.findall(r"POSITION\s+TOTAL-FORCE \(eV/Angst\)\s*\n[ \t]*-{3,}[ \t]*\n(.*?)(?=\n[ \t]*-{3,}[ \t]*\n)", tail, re.S)
    if not blocks:
        raise ValueError("No final OUTCAR positions/forces")
    return np.array([[float(x) for x in line.split()] for line in blocks[-1].strip().splitlines()])


def normalize_setting(key, value):
    value = value.upper().strip(".").replace("TRUE", "T").replace("FALSE", "F")
    # VASP 5 prints the executed METAGGA name in a five-character field.
    return "R2SCAN" if key == "METAGGA" and value == "R2SCA" else value


def audit(directory, surface, molecule):
    row = {"surface": surface, "molecule": molecule, "functional": FUNCTIONALS[directory.parent.name],
           "source_directory": str(directory), "accepted": False}
    try:
        paths = {name: directory / name for name in ("CONTCAR", "OUTCAR", "INCAR", "POSCAR")}
        for path in paths.values():
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"Missing/empty {path.name}")
        settings = incar_values(paths["INCAR"].read_text())
        row["incar"] = settings
        if settings.get("SYSTEM") != f"{surface}_{molecule}":
            raise ValueError("INCAR SYSTEM does not match exact target")
        if int(settings.get("NSW", "0")) <= 0 or int(settings.get("IBRION", "-1")) not in (1, 2, 3):
            raise ValueError("INCAR does not enable ionic relaxation")
        functional = directory.parent.name
        if functional == "beef_vdw":
            valid_func = settings.get("GGA", "").upper() == "BF" and settings.get("LUSE_VDW", "").upper() in (".TRUE.", "T")
        elif functional == "r2scan":
            valid_func = settings.get("METAGGA", "").upper() == "R2SCAN"
        else:
            valid_func = settings.get("GGA", "PE").upper() == "PE" and not settings.get("METAGGA")
            valid_func = valid_func and (settings.get("IVDW", "0") in ("11", "12") if functional == "PBE_D3" else settings.get("IVDW", "0") == "0")
        if not valid_func:
            raise ValueError("Functional directory disagrees with INCAR")
        output = read_output(paths["OUTCAR"])
        row["convergence"] = output["convergence"]
        if output["convergence"] != "converged":
            raise ValueError(output["convergence"])
        atoms = read(paths["CONTCAR"], format="vasp")
        initial = read(paths["POSCAR"], format="vasp")
        counts = Counter(atoms.get_chemical_symbols())
        metal = re.match(r"[A-Z][a-z]?", surface).group()
        if counts != Counter(initial.get_chemical_symbols()) or dict(counts) != output["composition"]:
            raise ValueError("CONTCAR/POSCAR/OUTCAR atom counts disagree")
        adsorbate = dict(counts)
        row["metal_count"] = adsorbate.pop(metal, 0)
        if row["metal_count"] <= 0 or adsorbate != FORMULAS[molecule]:
            raise ValueError("Wrong metal or molecular composition")
        row["composition"] = dict(counts)
        if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell).all():
            raise ValueError("Non-finite CONTCAR geometry")
        with paths["OUTCAR"].open("rb") as source:
            head = source.read(100000).decode(errors="replace")
            source.seek(max(0, paths["OUTCAR"].stat().st_size - 1000000))
            tail = source.read().decode(errors="replace")
        # Check the executed settings as well as the input file, including stale inputs.
        for key in ("NSW", "IBRION", "EDIFFG", "EDIFF"):
            executed = re.findall(r"\b" + key + r"\s*=\s*([-+\d.Ee]+)", head)
            if not executed or not np.isclose(float(executed[-1]), float(settings[key]), rtol=1e-6, atol=1e-12):
                raise ValueError(f"Executed {key} disagrees with INCAR")
        for key in ("GGA", "METAGGA", "IVDW", "LUSE_VDW"):
            if key not in settings:
                continue
            executed = re.findall(r"\b" + key + r"\s*=\s*([^\s;]+)", head)
            if not executed or normalize_setting(key, executed[-1]) != normalize_setting(key, settings[key]):
                raise ValueError(f"Executed {key} disagrees with INCAR")
        data = final_positions_forces(tail)
        if data.shape != (len(atoms), 6) or not np.isfinite(data).all():
            raise ValueError("Invalid final OUTCAR positions/forces")
        lattices = re.findall(r"direct lattice vectors\s+reciprocal lattice vectors\s*\n((?:[^\n]+\n){3})", head + tail)
        lattice = np.array([[float(x) for x in line.split()[:3]] for line in lattices[-1].strip().splitlines()])
        if not np.allclose(lattice, atoms.cell, atol=2e-5, rtol=0):
            raise ValueError("CONTCAR cell differs from OUTCAR")
        _, displacement = find_mic(atoms.positions - data[:, :3], atoms.cell, pbc=True)
        row["max_position_difference_A"] = float(max(displacement))
        if max(displacement) > 2e-5:
            raise ValueError("CONTCAR differs from final OUTCAR positions")
        forces = data[:, 3:].copy()
        for constraint in atoms.constraints:
            constraint.adjust_forces(atoms, forces)
        row["max_unconstrained_force_eV_A"] = float(np.linalg.norm(forces, axis=1).max())
        row["ediffg_eV_A"] = float(settings["EDIFFG"])
        if row["ediffg_eV_A"] < 0 and row["max_unconstrained_force_eV_A"] > abs(row["ediffg_eV_A"]) + 1e-5:
            raise ValueError("Final unconstrained force exceeds EDIFFG")
        row["ionic_steps"] = len(re.findall(r"Iteration\s+(\d+)\(\s*1\)", paths["OUTCAR"].read_text(errors="replace")))
        row["sources"] = {name: {"path": str(path), "sha256": digest(path)} for name, path in paths.items()}
        row["accepted"] = True
    except (ValueError, OSError, KeyError, IndexError) as error:
        row["reason"] = str(error)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--audit-output", type=Path)
    args = parser.parse_args()
    page_path = ROOT / "dft_comparison.html"
    page = page_path.read_text()
    soup = BeautifulSoup(page, "html.parser")
    missing = {(card["data-surf"], card["data-mol"]) for card in soup.select(".g")
               if card.select_one(".pair > div:nth-child(2) .ph")}
    if not missing:
        print("No missing DFT placeholders; no changes.")
        return
    candidates = defaultdict(list)
    visited = 0
    for parent, directories, _ in os.walk(args.source_root.resolve()):
        if "fully_relaxed" not in directories:
            continue
        directory = Path(parent) / "fully_relaxed"
        visited += 1
        if directory.parent.name not in FUNCTIONALS:
            continue
        # Find an exact target ancestor, skipping singlepoint and carbon-group folders.
        for ancestor in directory.parents:
            pair = tuple(ancestor.name.split("_", 1))
            if pair in missing:
                candidates[pair].append(directory)
                break
    rows, selected = [], []
    for surface, molecule in sorted(missing):
        checked = [audit(path, surface, molecule) for path in sorted(candidates[(surface, molecule)])]
        valid = [row for row in checked if row["accepted"]]
        # Honor the supplied Ag111 formate example; elsewhere use the existing renderer priority.
        priority = (["BEEF-vdW", "PBE", "PBE+D3", "r2SCAN"] if (surface, molecule) == ("Ag111", "formate")
                    else list(FUNCTIONALS.values()))
        if valid:
            chosen = min(valid, key=lambda row: (priority.index(row["functional"]), len(Path(row["source_directory"]).parts), row["source_directory"]))
            chosen["selected"] = True
            chosen["image"] = f"dftcmp/png/{surface}_{molecule}_dft.png"
            chosen["available_converged_functionals"] = sorted({row["functional"] for row in valid})
            selected.append(chosen)
        rows.extend(checked)
        print(f"{surface}_{molecule}: {len(valid)}/{len(checked)} verified; " + (chosen["functional"] if valid else "no image"), flush=True)
    report = {"source_root": str(args.source_root.resolve()), "fully_relaxed_directories_scanned": visited,
              "rotation": ROTATION, "show_unit_cell": 2, "ase_version": ase.__version__,
              "selection_policy": "Ag111 formate: supplied BEEF-vdW example first. Otherwise PBE, PBE+D3, r2SCAN, BEEF-vdW. Ties: shallowest path, then lexical path. Only verified completed relaxations qualify; energies are not used.",
              "missing_targets": [list(pair) for pair in sorted(missing)], "candidates": rows}
    if args.audit_output:
        args.audit_output.write_text(json.dumps(report, indent=2) + "\n")
    if not args.apply:
        return
    replacements = {}
    for chosen in selected:
        surface, molecule = chosen["surface"], chosen["molecule"]
        image_path = ROOT / chosen["image"]
        if image_path.exists():
            raise ValueError(f"Refusing to overwrite existing image: {image_path}")
        atoms = read(chosen["sources"]["CONTCAR"]["path"], format="vasp")
        write(image_path, atoms, rotation=ROTATION, show_unit_cell=2)
        chosen["image_sha256"] = digest(image_path)
        label = html.escape(chosen["functional"])
        source = html.escape(chosen["sources"]["CONTCAR"]["path"], quote=True)
        replacements[(surface, molecule)] = (f'<div><img src="{chosen["image"]}" alt="{surface} {molecule} DFT {label} relaxed" '
            f'title="Source: {source}"><div class="cap">DFT (relaxed, {label}) &middot; '
            f'<a href="dft_structure_sources.json" title="Source paths, functional selection and convergence checks">source audit</a></div></div>')
    placeholder = '<div><div class="ph">DFT structure<br>pending cluster<br>extraction</div><div class="cap">DFT (relaxed)</div></div>'
    def update_card(match):
        key = tuple(html.unescape(x) for x in match.group(2, 3))
        block = match.group(0)
        if key in replacements:
            assert block.count(placeholder) == 1
            return block.replace(placeholder, replacements[key])
        return block
    updated = re.sub(r'(<div class="g" data-surf="([^"]+)" data-mol="([^"]+)">).*?(?=\n   <div class="g"|\n</div>)', update_card, page, flags=re.S)
    count = len(soup.select('.pair > div:nth-child(2) img')) + len(selected)
    updated = re.sub(r'DFT IMAGES \d+/415', f'DFT IMAGES {count}/415', updated)
    updated = updated.replace('Structure images are from the earlier geometry extraction;',
        'Existing structure images are from the earlier geometry extraction; added Kestrel images carry their own functional and source audit. Geometry metrics retain the earlier extraction, and an image does not imply a new adsorption-energy result;')
    # Tables must remain byte-for-byte identical, including all values and provenance tooltips.
    assert re.findall(r'<table\b.*?</table>', updated, re.S) == re.findall(r'<table\b.*?</table>', page, re.S)
    parsed = BeautifulSoup(updated, "html.parser")
    for chosen in selected:
        assert parsed.select_one(f'img[src="{chosen["image"]}"]')
    (ROOT / "dft_structure_sources.json").write_text(json.dumps(report, indent=2) + "\n")
    page_path.write_text(updated)
    print(f"Added {len(selected)} images; all energy tables preserved byte-for-byte.")


if __name__ == "__main__":
    main()
