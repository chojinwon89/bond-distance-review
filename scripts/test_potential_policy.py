import csv,unittest
from pathlib import Path
from bs4 import BeautifulSoup
from potential_matching import potential_match

class PotentialPolicyTests(unittest.TestCase):
    def test_identity_requires_variant_and_date(self):
        c={'potentials':{'Pd':'PAW_PBE Pd 04Jan2005'}}
        for identity,expected in [('PAW_PBE Pd 04Jan2005',True),('PAW_PBE Pd_pv 28Jan2005',False),('PAW_PBE Pd 01Jan2020',False),(None,False)]:
            self.assertEqual(potential_match(c,{'composition':{'Pd':36},'potentials':{'Pd':identity}}),expected)

    def test_every_displayed_energy_is_a_selected_matching_candidate(self):
        root=Path(__file__).resolve().parents[1]
        if not (root/'dft_potential_mismatches.csv').exists():self.skipTest('Strict audit not generated')
        from update_kestrel_singlepoint import FUNCTIONALS
        from molecule_names import canonical
        with (root/'dft_cluster_selection.csv').open() as f:rows={(r['surface'],r['molecule'],r['functional'],r['mode']):r for r in csv.DictReader(f)}
        page=BeautifulSoup((root/'dft_comparison.html').read_text(),'html.parser')
        for card in page.select('.g'):
            for tr in card.select('table tr'):
                cells=tr.select('td')
                if not cells:continue
                for mode,index in [('relaxed',1),('SPE',4)]:
                    try:value=float(cells[index].get_text())
                    except ValueError:continue
                    key=(card['data-surf'],canonical(card['data-mol']),FUNCTIONALS[cells[0].get_text()],mode)
                    selected=rows[key];self.assertEqual(selected['application_status'],'selected candidate applied',key)
                    self.assertAlmostEqual(value,float(selected['E_ads']),delta=0.00051,msg=str(key))

if __name__=='__main__':unittest.main()
