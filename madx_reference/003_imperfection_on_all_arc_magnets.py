import numpy as np
from scipy.special import factorial

import xtrack as xt

line = xt.load('./hllhc_b1_v19_round_imo300_no_errors.json')

# Load WISE error table
tt_raw = xt.Table.from_tfs('madx/errors/LHC/wise/collision_errors-emfqcs-6.tfs')
tt_err_data = tt_raw._data.copy()
tt_err_data['name'] = np.array([nn.lower() for nn in tt_err_data['name']])
tt_err = xt.Table(data=tt_err_data, col_names=tt_raw._col_names)

# dicts with LHC geography
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


# Isolate magnets with two apertures (end with .b1, .b2, .v1, .v2)is_two_aper = np.array([nn.endswith('.b1') or nn.endswith('.b2') or nn.endswith('.v1') or nn.endswith('.v2') for nn in line.element_names])
tt_err_two_aper = tt_err.rows[r'.*\.b1|.*\.b2|.*\.v1|.*\.v2']

# Use name with beam instead of name with aper
name_with_aper = tt_err_two_aper['name']
name_with_beam = []
for nn in name_with_aper:
    side_aper = nn[-5:]
    side_beam = SIDE_APER_TO_SIDE_BEAM[side_aper]
    name_with_beam.append(nn[:-5] + side_beam)
name_with_beam = np.array(name_with_beam)

tt_err_two_aper['name_with_aper'] = name_with_aper
tt_err_two_aper['name_with_beam'] = name_with_beam
tt_err_two_aper['name'] = name_with_beam

# Attach reference order (0 if mb, 1 if mq, fail otherwise)
ref_order = []
for nn in tt_err_two_aper['name']:
    if nn.startswith('mb'):
        ref_order.append(0)
    elif nn.startswith('mq'):
        ref_order.append(1)
    else:
        raise ValueError(f"Unexpected magnet name: {nn}")
tt_err_two_aper['ref_order'] = np.array(ref_order)

# From WISE units to knl_rel and ksl_rel
max_order = 15
ref_radius = 0.017

knl_rel = np.zeros((len(tt_err_two_aper), max_order))
ksl_rel = np.zeros((len(tt_err_two_aper), max_order))

for ii in range(0, max_order):

    aa = tt_err_two_aper[f'a{ii+1}']
    bb = tt_err_two_aper[f'b{ii+1}']
    ref_order = tt_err_two_aper['ref_order']

    dknlr_mad = 1e-4 * bb * (-1) ** (ii    )
    dkslr_mad = 1e-4 * aa * (-1) ** (ii + 1)

    kknn_rel = dknlr_mad * ref_radius**(ref_order - (ii)) * factorial(ii) / factorial(ref_order)
    kkss_ref = dkslr_mad * ref_radius**(ref_order - (ii)) * factorial(ii) / factorial(ref_order)

    knl_rel[:, ii] = kknn_rel
    ksl_rel[:, ii] = kkss_ref

tt_err_two_aper['knl_rel'] = knl_rel
tt_err_two_aper['ksl_rel'] = ksl_rel