import numpy as np

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
tt_two_aper = tt_err.rows[r'.*\.b1|.*\.b2|.*\.v1|.*\.v2']

# Attach name with beam instead of aperture to error table

name_with_aper = tt_two_aper['name']
name_with_beam = np.array([SIDE_APER_TO_SIDE_BEAM.get(nn, nn) for nn in name_with_aper])
tt_two_aper['name_with_aper'] = name_with_aper
tt_two_aper['name_with_beam'] = name_with_beam
tt_two_aper['name'] = name_with_beam


tt = line.get_table()
name_no_slice = []
slice_index = []
for nn in tt.name:
    if '..' in nn and not nn.endswith('..r'): # ignore the "r" correctors ?!
        nn_no_slice = nn.split('..')[0]
        name_no_slice.append(nn_no_slice)
        slice_index.append(int(nn.split('..')[1]))
    else:
        name_no_slice.append(nn)
        slice_index.append(None)
