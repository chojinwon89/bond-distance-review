"""Require identical recorded PAW TITEL identities for every reference species."""
def potential_match(complex_result, reference):
    cp=complex_result.get('potentials',{});rp=reference.get('potentials',{})
    return bool(reference.get('composition')) and all(
        cp.get(element) and cp.get(element)==rp.get(element)
        for element in reference['composition'])
