import xtrack as xt
import numpy as np
from scipy.special import factorial

from lhc_geography import SIDE_APER_TO_SIDE_BEAM

def load_wise_table_arc_magnets(fname, min_order=2, max_order=15, ref_radius=0.017):
    # Load WISE error table
    tt_raw = xt.Table.from_tfs(fname)
    tt_err_data = tt_raw._data.copy()
    tt_err_data['name'] = np.array([nn.lower() for nn in tt_err_data['name']])
    tt_err = xt.Table(data=tt_err_data, col_names=tt_raw._col_names)

    # Isolate magnets with two apertures (end with .b1, .b2, .v1, .v2)
    tt_err_two_aper = tt_err.rows[r'.*\.b1|.*\.b2|.*\.v1|.*\.v2']

    # I want to keep only the arcs (skip cells 1-7)
    for icell in range(1, 8):
        for ip in [1, 2, 3, 4, 5, 6, 7, 8]:
            tt_err_two_aper = tt_err_two_aper.rows.match_not(
                f'.*\\.{icell}r{ip}.*|.*\\.{icell}l{ip}.*'
                f'|.*\\.a{icell}r{ip}.*|.*\\.a{icell}l{ip}.*'
                f'|.*\\.b{icell}r{ip}.*|.*\\.b{icell}l{ip}.*'
            )

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
    ref_radius = 0.017

    knl_rel = np.zeros((len(tt_err_two_aper), max_order))
    ksl_rel = np.zeros((len(tt_err_two_aper), max_order))

    for ii in range(0, max_order):

        aa = tt_err_two_aper[f'a{ii+1}']
        bb = tt_err_two_aper[f'b{ii+1}']
        ref_order = tt_err_two_aper['ref_order']

        dknlr_mad = 1e-4 * bb * (-1) ** (ref_order + ii    )
        dkslr_mad = 1e-4 * aa * (-1) ** (ref_order + ii + 1)

        kknn_rel = dknlr_mad * ref_radius**(ref_order - (ii)) * factorial(ii) / factorial(ref_order)
        kkss_ref = dkslr_mad * ref_radius**(ref_order - (ii)) * factorial(ii) / factorial(ref_order)

        knl_rel[:, ii] = kknn_rel
        ksl_rel[:, ii] = kkss_ref

    tt_err_two_aper['knl_rel'] = knl_rel
    tt_err_two_aper['ksl_rel'] = ksl_rel

    multipole_errors = {}
    for nn in tt_err_two_aper.name:
        knl_rel = tt_err_two_aper['knl_rel', nn]
        ksl_rel = tt_err_two_aper['ksl_rel', nn]
        ref_order = tt_err_two_aper['ref_order', nn]
        multipole_errors[nn] = {'knl_rel': knl_rel, 'ksl_rel': ksl_rel, 'ref_order': ref_order}

    # Suppress multipoles of order < 2
    for nn in multipole_errors:
        multipole_errors[nn]['knl_rel'][:min_order] = 0
        multipole_errors[nn]['ksl_rel'][:min_order] = 0

    return multipole_errors