import gzip,hashlib,json,tempfile,unittest
from pathlib import Path
from component_store import read_store
from cluster_result_bundle import import_bundle
from test_shared_reference_energies import components

class BundleTests(unittest.TestCase):
    def make_bundle(self,root,corrupt=False):
        rows=components();c=rows[0];text='saved VASP coordinates\n'
        payload=dict(schema_version=1,cluster='perlmutter',components=rows,coordinates=[dict(cluster=c['cluster'],directory=c['calculation']['directory'],snapshot_id=c['snapshot_id'],sha256=hashlib.sha256(text.encode()).hexdigest(),vasp=text+('corrupt' if corrupt else ''))])
        path=root/'bundle.json.gz';path.write_bytes(gzip.compress(json.dumps(payload).encode()));return path

    def test_merge_is_idempotent_and_caches_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);path=self.make_bundle(root);import_bundle(path,root);first=(root/'dft_component_store.jsonl.gz').read_bytes();import_bundle(path,root)
            self.assertEqual(first,(root/'dft_component_store.jsonl.gz').read_bytes());self.assertEqual(len(read_store(root/'dft_component_store.jsonl.gz')),3)
            index=json.loads((root/'data/adsorption_results/coordinate_sources.json').read_text());self.assertEqual(len(index),1)
            self.assertTrue((root/next(iter(index.values()))['path']).exists())

    def test_invalid_coordinate_hash_prevents_store_mutation(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);path=self.make_bundle(root,True)
            with self.assertRaises(ValueError):import_bundle(path,root)
            self.assertFalse((root/'dft_component_store.jsonl.gz').exists())

if __name__=='__main__':unittest.main()
