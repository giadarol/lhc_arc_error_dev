import xtrack as xt

from load_wise import load_wise_table_arc_magnets

min_order = 2
max_order = 15

line = xt.load('./hllhc_b1_v19_round_imo300_no_errors.json')

multipole_errors = load_wise_table_arc_magnets(
    'madx/errors/LHC/wise/collision_errors-emfqcs-6.tfs',
    min_order=min_order, max_order=max_order)

# Apply errors in the line
for nn in line.element_names:
    if not hasattr(line[nn], 'knl'):
        continue  # skip non-multipoles

    print(f'Applying errors to {nn}               ', end='\r', flush=True)
    nn_err = nn.split('..')[0]  # remove ..1, ..2, etc.
    if nn_err in multipole_errors:
        line.extend_knl_ksl(order=max_order, element_names=[nn])
        for ii in range(min_order, max_order):
            kknn_rel = multipole_errors[nn_err]['knl_rel'][ii]
            kkss_rel = multipole_errors[nn_err]['ksl_rel'][ii]
            ref_order = int(multipole_errors[nn_err]['ref_order'])
            if ii == ref_order:
                raise ValueError(
                    f"Error of order {ii} for {nn} is relative to the reference multipole, which is not supported. "
                )
            line.ref[nn].knl[ii] += kknn_rel * line.ref[nn].knl[ref_order]
            line.ref[nn].ksl[ii] += kkss_rel * line.ref[nn].knl[ref_order]
            # line.get(nn).knl[ii] = kknn_rel * line.get(nn).knl[ref_order]
            # line.get(nn).ksl[ii] = kkss_rel * line.get(nn).knl[ref_order]

line.to_json('test_line_with_errors.json')