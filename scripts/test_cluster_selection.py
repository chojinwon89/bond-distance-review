import unittest
from cluster_selection import choose
from structure_identity import same_structure


def candidate(cluster,status,geometry=None):
    c=dict(directory='/'+cluster+'/DME/PBE',composition={'Ag':36,'C':2,'H':6,'O':1},cell=[[1,0,0],[0,1,0],[0,0,1]])
    if geometry:c['structure_identity']={'input':{'geometry_sha256':geometry},'relaxed':{'geometry_sha256':geometry}}
    return dict(status=status,complex_cluster=cluster,_complex={'calculation':c})


class ClusterPolicyTests(unittest.TestCase):
    def test_normal_perlmutter_wins(self):
        p=candidate('perlmutter','ok','same');k=candidate('kestrel','ok','same')
        self.assertIs(choose(p,k,'SPE')[0],p)

    def test_verified_kestrel_fallback(self):
        p=candidate('perlmutter','energy_review','same');k=candidate('kestrel','ok','same')
        selected,reason,fallback=choose(p,k,'SPE')
        self.assertIs(selected,k);self.assertTrue(fallback);self.assertIn('same structure verified',reason)

    def test_names_and_cells_do_not_prove_same_structure(self):
        for left,right in [(None,None),('one','two')]:
            p=candidate('perlmutter','energy_review',left);k=candidate('kestrel','ok',right)
            self.assertIs(choose(p,k,'SPE')[0],p)
            self.assertFalse(same_structure(p['_complex'],k['_complex'],'SPE')[0])

    def test_relaxed_match_requires_final_geometry(self):
        p=candidate('perlmutter','energy_review','same');k=candidate('kestrel','ok','same')
        del k['_complex']['calculation']['structure_identity']['relaxed']
        self.assertIs(choose(p,k,'relaxed')[0],p)
        self.assertIs(choose(p,k,'SPE')[0],k)

    def test_two_unusual_candidates_keep_perlmutter(self):
        p=candidate('perlmutter','energy_review','same');k=candidate('kestrel','energy_review','same')
        self.assertIs(choose(p,k,'SPE')[0],p)


if __name__=='__main__':unittest.main()
