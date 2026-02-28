import xtrack as xt
import numpy as np
from math import factorial

# TODO:
# - Remember to mask away low-order multipoles
# - yrot, s_rot where applicable
# - b2

line_ref = xt.load("hllhc_b1_v19_round_imo300_arcs.json")

APER_NAME = {
    'b1': {
        '12': 'v1',
        '23': 'v2',
        '34': 'v2',
        '45': 'v2',
        '56': 'v1',
        '67': 'v1',
        '78': 'v1',
        '81': 'v2',
    },
    'b2': {
        '12': 'v2',
        '23': 'v1',
        '34': 'v1',
        '45': 'v1',
        '56': 'v2',
        '67': 'v2',
        '78': 'v2',
        '81': 'v1',
    }
}

tt_raw = xt.Table.from_tfs('madx/errors/LHC/wise/collision_errors-emfqcs-6.tfs')
tt_err_data = tt_raw._data.copy()
tt_err_data['name'] = np.array([nn.lower() for nn in tt_err_data['name']])
tt_err = xt.Table(data=tt_err_data, col_names=tt_raw._col_names)

beam = 'b1'
sector = '56'
slice_name = 'mb.a25r5.b1..1'
magnet = slice_name.split('..')[0]
order_ref = 0

aper_name = '.'.join(magnet.split('.')[:-1]) + '.' + APER_NAME[beam][sector]

ref_r = 0.017
max_order = 15

dknlr_mad = []
dkslr_mad = []

for ii in range(0, max_order):
    aa = tt_err[f'a{ii+1}', aper_name]
    bb = tt_err[f'b{ii+1}', aper_name]

    dknlr_mad.append(1e-4 * bb)
    dkslr_mad.append(-1e-4 * aa)

dklnr_mad = np.array(dknlr_mad)
dksnr_mad = np.array(dkslr_mad)

multipole_obj = line_ref[slice_name]
k_ref = multipole_obj.knl[order_ref]

knl = []
ksl = []

for ii in range(0, max_order):
    kknn = dklnr_mad[ii] * k_ref * ref_r**(order_ref - (ii)) * factorial(ii) / factorial(order_ref)
    kkss = dksnr_mad[ii] * k_ref * ref_r**(order_ref - (ii)) * factorial(ii) / factorial(order_ref)
    knl.append(kknn)
    ksl.append(kkss)

knl = np.array(knl)
ksl = np.array(ksl)
