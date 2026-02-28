import numpy as np
from scipy.special import factorial

import xtrack as xt

ERROR_KNOBS = {
    "ON_A1s": 0, "ON_A1r": 0, "ON_B1s": 0, "ON_B1r": 0, 
    "ON_A2s": 0, "ON_A2r": 0, "ON_B2s": 0, "ON_B2r": 0, 
    "ON_A3s": 1, "ON_A3r": 1, "ON_B3s": 1, "ON_B3r": 1, 
    "ON_A4s": 1, "ON_A4r": 1, "ON_B4s": 1, "ON_B4r": 1, 
    "ON_A5s": 1, "ON_A5r": 1, "ON_B5s": 1, "ON_B5r": 1, 
    "ON_A6s": 1, "ON_A6r": 1, "ON_B6s": 1, "ON_B6r": 1, 
    "ON_A7s": 1, "ON_A7r": 1, "ON_B7s": 1, "ON_B7r": 1, 
    "ON_A8s": 1, "ON_A8r": 1, "ON_B8s": 1, "ON_B8r": 1, 
    "ON_A9s": 1, "ON_A9r": 1, "ON_B9s": 1, "ON_B9r": 1, 
    "ON_A10s": 1, "ON_A10r": 1, "ON_B10s": 1, "ON_B10r": 1, 
    "ON_A11s": 1, "ON_A11r": 1, "ON_B11s": 1, "ON_B11r": 1, 
    "ON_A12s": 1, "ON_A12r": 1, "ON_B12s": 1, "ON_B12r": 1, 
    "ON_A13s": 1, "ON_A13r": 1, "ON_B13s": 1, "ON_B13r": 1, 
    "ON_A14s": 1, "ON_A14r": 1, "ON_B14s": 1, "ON_B14r": 1, 
    "ON_A15s": 1, "ON_A15r": 1, "ON_B15s": 1, "ON_B15r": 1, 
    "ON_A16s": 1, "ON_A16r": 1, "ON_B16s": 1, "ON_B16r": 1, 
    "ON_A17s": 1, "ON_A17r": 1, "ON_B17s": 1, "ON_B17r": 1, 
    "ON_A18s": 1, "ON_A18r": 1, "ON_B18s": 1, "ON_B18r": 1, 
    "ON_A19s": 1, "ON_A19r": 1, "ON_B19s": 1, "ON_B19r": 1, 
    "ON_A20s": 1, "ON_A20r": 1, "ON_B20s": 1, "ON_B20r": 1, 
}

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

SECTORS = {
    'r1': '12',
    'l2': '12',
    'r2': '23',
    'l3': '23',
    'r3': '34',
    'l4': '34',
    'r4': '45',
    'l5': '45',
    'r5': '56',
    'l6': '56', 
    'r6': '67', 
    'l7': '67', 
    'r7': '78', 
    'l8': '78', 
    'r8': '81', 
    'l1': '81'
}

ref_r = 0.017
max_order = 15


beam = 1
seed = 6
version = 19
optics = 'round'
i_mo = 300.
out_path = './'


# Load clean line
line = xt.load(f"{out_path}/hllhc_b{beam}_v{version}_{optics}_imo{i_mo:.0f}_no_errors.json")

# Load line with errors to compare
line_ref = xt.load(f"{out_path}/hllhc_b{beam}_v{version}_{optics}_imo{i_mo:.0f}_arcs.json")

# Get energy and load correct error table
nrj = line.particle_ref.p0c[0]
if nrj < 5000e9:
    table_path = f"madx/errors/LHC/wise/injection_errors-emfqcs-{seed}.tfs"
else:
    table_path = f"madx/errors/LHC/wise/collision_errors-emfqcs-{seed}.tfs"
tt_raw = xt.Table.from_tfs(table_path)
tt_err_data = tt_raw._data.copy()
tt_err_data['name'] = np.array([nn.lower() for nn in tt_err_data['name']])
tt_err = xt.Table(data=tt_err_data, col_names=tt_raw._col_names)

# Loop over all bends and install errors
elements, names = line.get_elements_of_type(xt.Multipole) # line must be thin
mask = np.strings.startswith(names, "mb.") # bend names start with MB.
elements = [e for e, m in zip(elements, mask) if m]
names = [n for n, m in zip(names, mask) if m]
line.env.extend_knl_ksl(max_order-1, names)
for name in names:
    # print(f"\nInstalling field errors for {name}")
    nn = name.split("..")[0] # elements are sliced
    sector = SECTORS[nn.split('.')[1][-2:]]
    aper_name = '.'.join(nn.split('.')[:-1]) + '.' + APER_NAME[f"b{beam}"][sector]

    order = 0
    k_ref = line[name].knl[0] / line[name].length

    dknlr_mad = []
    dkslr_mad = []

    for ii in range(0, max_order):
        aa = tt_err[f'a{ii+1}', aper_name]
        bb = tt_err[f'b{ii+1}', aper_name]

        dknlr_mad.append((-1)**ii * bb * ERROR_KNOBS[f"ON_B{ii+1}s"] * 1e-4)
        dkslr_mad.append((-1)**(ii+1) * aa * ERROR_KNOBS[f"ON_A{ii+1}s"] * 1e-4)
    
    dknlr_mad = np.asarray(dknlr_mad) * line[name].length
    dkslr_mad = np.asarray(dkslr_mad) * line[name].length

    # Make it absolute error
    dknl = dknlr_mad * k_ref * ref_r**(order-np.arange(len(dknlr_mad))) * factorial(np.arange(len(dknlr_mad))) / factorial(order) / (2*(len(dknlr_mad)-np.arange(len(dknlr_mad)))+1)
    dksl = dkslr_mad * k_ref * ref_r**(order-np.arange(len(dkslr_mad))) * factorial(np.arange(len(dkslr_mad))) / factorial(order) / (2*(len(dkslr_mad)-np.arange(len(dkslr_mad))))

    for jj in range(len(dknl)):
        line.element_refs[name].knl[jj] = line.element_refs[name].knl[jj] + dknl[jj]
        line.element_refs[name].ksl[jj] = line.element_refs[name].ksl[jj] + dksl[jj]

    if np.any(np.where(np.abs((line[name].knl-line_ref[name].knl)/line_ref[name].knl*100) > 1e-10)):
        print(f"Something wrong with knl for {name}.")
        print(line[name].knl)
        print(line_ref[name].knl)
    if np.any(np.where(np.abs((line[name].ksl-line_ref[name].ksl)/line_ref[name].ksl*100) > 1e-10)):
        print(f"Something wrong with ksl for {name}.")
        print(line[name].ksl)
        print(line_ref[name].ksl)

# print(f"\n mb.a25r5.b1..1")
# print(line["mb.a25r5.b1..1"].knl)
# print(line_ref["mb.a25r5.b1..1"].knl)
# print((line["mb.a25r5.b1..1"].knl-line_ref["mb.a25r5.b1..1"].knl)/line_ref["mb.a25r5.b1..1"].knl*100)
# print(line["mb.a25r5.b1..1"].knl/line_ref["mb.a25r5.b1..1"].knl)
# print(line["mb.a25r5.b1..1"].ksl)
# print(line_ref["mb.a25r5.b1..1"].ksl)
# print((line["mb.a25r5.b1..1"].ksl-line_ref["mb.a25r5.b1..1"].ksl)/line_ref["mb.a25r5.b1..1"].ksl*100)
# print(line["mb.a25r5.b1..1"].ksl/line_ref["mb.a25r5.b1..1"].ksl)