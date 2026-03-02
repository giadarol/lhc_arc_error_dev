
B1_INSIDE = {'b1': 'v2', 'b2': 'v1',
             'v1': 'b2', 'v2': 'b1',
             'internal_beam': 'b1', 'external_beam': 'b2',
             'internal_aper': 'v2', 'external_aper': 'v1'}
B2_INSIDE = {'b1': 'v1', 'b2': 'v2',
             'v1': 'b1', 'v2': 'b2',
             'internal_beam': 'b2', 'external_beam': 'b1',
             'internal_aper': 'v2', 'external_aper': 'v1'}
B1_OUTSIDE = B2_INSIDE
B2_OUTSIDE = B1_INSIDE

BEAM_MAPPING_PER_SIDE = {
    'r1': B1_OUTSIDE, 'l2': B1_OUTSIDE,
    'r2': B1_INSIDE,  'l3': B1_INSIDE,
    'r3': B1_INSIDE,  'l4': B1_INSIDE,
    'r4': B1_INSIDE,  'l5': B1_INSIDE,
    'r5': B1_OUTSIDE, 'l6': B1_OUTSIDE,
    'r6': B1_OUTSIDE, 'l7': B1_OUTSIDE,
    'r7': B1_OUTSIDE, 'l8': B1_OUTSIDE,
    'r8': B1_INSIDE,  'l1': B1_INSIDE,
}

SIDE_APER_TO_SIDE_BEAM = {}
for side in BEAM_MAPPING_PER_SIDE:
    for aper in ['v1', 'v2']:
        SIDE_APER_TO_SIDE_BEAM[f"{side}.{aper}"] = f"{side}.{BEAM_MAPPING_PER_SIDE[side][aper]}"
