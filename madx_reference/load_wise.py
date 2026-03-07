import xtrack as xt
import numpy as np
from scipy.special import factorial

from lhc_geography import SIDE_APER_TO_SIDE_BEAM

def load_wise_table_arc_magnets(fname_err_table, fname_rotations, min_order=2, max_order=15, ref_radius=0.017):
    # Load WISE error table
    tt_raw = xt.Table.from_tfs(fname_err_table)
    tt_err_data = tt_raw._data.copy()
    tt_err_data['name'] = np.array([nn.lower() for nn in tt_err_data['name']])
    tt_err = xt.Table(data=tt_err_data, col_names=tt_raw._col_names)

    # Load rotation table
    tt_raw_rot = xt.Table.from_tfs(fname_rotations)
    tt_rot_data = tt_raw_rot._data.copy()
    tt_rot_data['name'] = np.array([nn.lower() for nn in tt_rot_data['name']])
    tt_rot = xt.Table(data=tt_rot_data, col_names=tt_raw_rot._col_names)
    rot={} # We turn it into a dict for easier access later.
    for nn in tt_rot['name']:
        rot[nn] = {'yrot': tt_rot['yrota', nn], 'srot': tt_rot['srota', nn], 'inout': tt_rot['inout', nn]}

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
    yrotfactor = []
    for nn in name_with_aper:
        # Check if in rot table
        nn_no_aper = nn.replace('.v1', '').replace('.v2', '')
        yyff = 1
        if nn_no_aper in rot: # aperture-beam mapping explicitly given in rot table
            inout = rot[nn_no_aper]['inout']
            if inout == 0:
                pass  # single bore, do nothing
            elif inout == 1:
                nn_with_beam = nn.replace('.v1', '.b1').replace('.v2', '.b2')
            elif inout == 2:
                nn_with_beam = nn.replace('.v1', '.b2').replace('.v2', '.b1')
            else:
                raise ValueError(f"Unexpected inout value {inout} for magnet {nn_no_aper}")
            yrot = rot[nn_no_aper]['yrot']
            assert yrot in [0, 180], f"Unexpected yrot value {yrot} for magnet {nn_no_aper}"
            if yrot == 180:
                yyff = -1
        else: # use standard LHC geography (v1 is the external beam)
            side_aper = nn[-5:]
            side_beam = SIDE_APER_TO_SIDE_BEAM[side_aper]
            nn_with_beam = nn[:-5] + side_beam
        name_with_beam.append(nn_with_beam)
        yrotfactor.append(yyff)
    name_with_beam = np.array(name_with_beam)

    tt_err_two_aper['name_with_aper'] = name_with_aper
    tt_err_two_aper['name_with_beam'] = name_with_beam
    tt_err_two_aper['name'] = name_with_beam
    tt_err_two_aper['yrotfactor'] = np.array(yrotfactor)

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

    ref_order = tt_err_two_aper['ref_order']
    yrotfactor = tt_err_two_aper['yrotfactor']
    for ii in range(0, max_order):

        aa = tt_err_two_aper[f'a{ii+1}']
        bb = tt_err_two_aper[f'b{ii+1}']

        # From magnet measurement convention to MADX convention
        dknlr_mad = 1e-4 * bb * (-1) ** (ref_order + ii    ) * yrotfactor ** ii
        dkslr_mad = 1e-4 * aa * (-1) ** (ref_order + ii + 1) * yrotfactor ** (ii + 1)

        # From MADX convention to knl
        kknn_rel = dknlr_mad * ref_radius**(ref_order - (ii)) * factorial(ii) / factorial(ref_order)
        kkss_rel = dkslr_mad * ref_radius**(ref_order - (ii)) * factorial(ii) / factorial(ref_order)

        knl_rel[:, ii] = kknn_rel
        ksl_rel[:, ii] = kkss_rel

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

    return multipole_errors, tt_err_two_aper