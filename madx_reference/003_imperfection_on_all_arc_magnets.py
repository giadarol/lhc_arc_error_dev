import xtrack as xt

from load_wise import load_wise_table_arc_magnets

# TODO:
# Handle skew
# Handle sign change for B2

# The rotations table is an invalid tfs because it has a space in a name
# We need to patch it (replace "not found" with "not_found")
fname_rotations = 'madx/errors/LHC/rotations_Q2_integral.tab'
with open(fname_rotations, 'r') as f:
    content = f.read()
content = content.replace('not found', 'not_found')
with open(fname_rotations + '_patched', 'w') as f:
    f.write(content)

min_order = 2
max_order = 15

line = xt.load('./hllhc_b1_v19_round_imo300_no_errors.json')

multipole_errors, tt_err_arc = load_wise_table_arc_magnets(
    fname_err_table='madx/errors/LHC/wise/collision_errors-emfqcs-6.tfs',
    fname_rotations=fname_rotations + '_patched',
    min_order=min_order, max_order=max_order)

element_names = line.element_names
error_knob_name = 'on_error_arc'
append_order_to_knob_name = True


env = line.env

if append_order_to_knob_name:
    for ii in range(min_order, max_order):
        env[f'{error_knob_name}_k{ii}'] = 1
        env[f'{error_knob_name}_k{ii}s'] = 1
else:
    env[error_knob_name] = 1

# Apply errors in the line
for nn in element_names:
    if not hasattr(line[nn], 'knl'):
        continue  # skip non-multipoles

    print(f'Applying errors to {nn}               ', end='\r', flush=True)
    nn_err = nn.split('..')[0]  # remove ..1, ..2, etc.
    if nn_err in multipole_errors:
        line.extend_knl_rel_ksl_rel(order=max_order, element_names=[nn])
        for ii in range(min_order, max_order):
            kknn_rel = multipole_errors[nn_err]['knl_rel'][ii]
            kkss_rel = multipole_errors[nn_err]['ksl_rel'][ii]
            ref_order = int(multipole_errors[nn_err]['ref_order'])
            if ii == ref_order:
                raise ValueError(
                    f"Error of order {ii} for {nn} is relative to the reference multipole, which is not supported. "
                )

            if append_order_to_knob_name:
                knob_kn_name = f'{error_knob_name}_k{ii}'
                knob_ks_name = f'{error_knob_name}_k{ii}s'
            else:
                knob_kn_name = error_knob_name
                knob_ks_name = error_knob_name

            # Using knl_rel and ksl_rel
            line[nn].main_order = ref_order
            line.ref[nn].knl_rel[ii] = kknn_rel * env.ref[knob_kn_name]
            line.ref[nn].ksl_rel[ii] = kkss_rel * env.ref[knob_ks_name]

line.to_json('test_line_with_errors.json')