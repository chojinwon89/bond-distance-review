"""Explicit chemical-name equivalences; never infer identity from a gross formula."""
SUBSCRIPTS=str.maketrans('₀₁₂₃₄₅₆₇₈₉','0123456789')
GROUPS={
    'DME':('CH3OCH3','dimethyl_ether','dimethyl ether','dme'),
    'methanol':('CH3OH','methyl_alcohol'),
    'ethanol':('CH3CH2OH','C2H5OH','ethyl_alcohol'),
    'acetaldehyde':('CH3CHO','ethanal'),
    'acetic_acid':('CH3COOH','ethanoic_acid','acetic acid'),
    'formaldehyde':('H2CO','CH2O','methanal'),
    'formic_acid':('HCOOH','methanoic_acid','formic acid'),
    'formate':('HCOO',),
    'methoxy':('CH3O',),
    'CH2OH':('hydroxymethyl',),
    'CH3':('methyl',),
    'methane':('CH4',),
    'ethane':('C2H6',),
    'ethene':('C2H4','ethylene'),
    'acetylene':('C2H2','ethyne'),
    'propane':('C3H8',),
    'propene':('propylene','CH3CH=CH2'),
    'propanol':('1-propanol','n-propanol','1_propanol','n_propanol','CH3CH2CH2OH'),
    'isopropanol':('2-propanol','2_propanol','isopropyl_alcohol','CH3CHOHCH3','(CH3)2CHOH'),
    'glycerol':('glycerin','HOCH2CHOHCH2OH'),
    'H2O':('water',),'NH3':('ammonia',),'CO':('carbon_monoxide',),
    'CO2':('carbon_dioxide',),'H2':('hydrogen',),'N2':('nitrogen',),
    'O2':('oxygen',),'H2S':('hydrogen_sulfide',),'SO2':('sulfur_dioxide',),
    'NO':('nitric_oxide',),
}
FORMULAS={key:values[0] for key,values in GROUPS.items() if key not in ('CH3','CH2OH','H2O','NH3','CO','CO2','H2','N2','O2','H2S','SO2','NO')}
FORMULAS.update(CH3='CH3',CH2OH='CH2OH',H2O='H2O',NH3='NH3',CO='CO',CO2='CO2',H2='H2',N2='N2',O2='O2',H2S='H2S',SO2='SO2',NO='NO',
                propene='CH3CH=CH2',propanol='CH3CH2CH2OH',isopropanol='CH3CHOHCH3',glycerol='HOCH2CHOHCH2OH')
ALIASES={}
for key,values in GROUPS.items():
    for name in (key,*values):
        if name in ALIASES and ALIASES[name]!=key:raise ValueError('Ambiguous alias: '+name)
        ALIASES[name]=key


def canonical(name):
    normalized=name.strip().translate(SUBSCRIPTS)
    return ALIASES.get(normalized,normalized)


def equivalent_names(name, available, preferred=None):
    """Prefer original spelling, then the prior reference spelling; never by energy."""
    return sorted((n for n in available if canonical(n)==canonical(name)),
                  key=lambda n:(n!=name,n!=preferred,n!=canonical(name),n))
