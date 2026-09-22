import tempfile
import unittest
from pathlib import Path
from publication_store import capture, read, REL


def page(value, attrs=''):
    return f'<div class="g" data-surf="Ag111" data-mol="CH3OCH3"><table><tr><td>PBE</td><td {attrs}>{value}</td><td>-0.20</td><td>0</td><td>—</td><td>—</td><td>—</td></tr></table></div>'


class PublicationStoreTests(unittest.TestCase):
    def test_failed_refresh_preserves_energy_and_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            capture(root,page('-0.123'),'git:original')
            selected=capture(root,page('refs unverified'),'local-before-selection')
            self.assertEqual(selected[('Ag111','DME','pbe','relaxed')]['energy_eV'],'-0.123')
            self.assertEqual(len(read(root)),1)

    def test_verified_update_retains_both_observations_idempotently(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            capture(root,page('-25'),'git:original')
            capture(root,page('-0.35'),'validated-selection')
            before=(root/REL).read_bytes()
            selected=capture(root,page('-0.35'),'validated-selection')
            self.assertEqual(before,(root/REL).read_bytes())
            self.assertEqual(len(read(root)),2)
            self.assertEqual(next(iter(selected.values()))['energy_eV'],'-0.35')
            selected=capture(root,page('-99'),'local-before-selection')
            self.assertEqual(next(iter(selected.values()))['energy_eV'],'-0.35')

    def test_historical_render_does_not_create_new_observation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            capture(root,page('-25'),'git:original')
            capture(root,page('-25','class="energy-review historical-energy"'),'validated-selection')
            self.assertEqual(len(read(root)),1)

    def test_corruption_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            capture(root,page('-25'),'git:original')
            p=root/REL;p.write_text(p.read_text().replace('"energy_eV": "-25"','"energy_eV": "-26"'))
            with self.assertRaises(ValueError):read(root)
