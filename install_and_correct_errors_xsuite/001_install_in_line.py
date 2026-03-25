import xtrack as xt
import xmask as xm
import xmask.lhc as xmlhc

# Load multipolar errors from json files
multipole_errors = xt.json.load('multipole_errors_arc.json')

# Get a line model
line = xt.load('hllhc_b1_v19_round_imo300_no_errors.json')
lhc = line.env
lhc['b1'] = line

# Apply errors in lines
min_order = 1
max_order = 15

xm.set_multipole_errors_in_line(line, multipole_errors,
                        min_order=min_order, max_order=max_order,
                        error_knob_name='on_error_arc',
                        append_order_to_knob_name=True)

# Switch off errors of order 0 and 1
lhc['on_error_arc_k0'] = 0
lhc['on_error_arc_k0s'] = 0
lhc['on_error_arc_k1'] = 0
lhc['on_error_arc_k1s'] = 0

line.vars.get_table().rows['on_error_.*'].show()

# Status of error knobs
tt_err_knobs = lhc.vars.get_table().rows[r'on_error_.*']
print("Error knobs in the environment:")
tt_err_knobs.show()

# Errors off to get reference twiss
lhc.set(tt_err_knobs.name, 0)
tw_b1 = lhc['b1'].twiss4d(reverse=False) # Reference twiss
# errors back on
for nn in tt_err_knobs.name:
    lhc[nn] = tt_err_knobs['value', nn]

# Spool piece correctors (MCS, MC0, MCD)
xmlhc.set_arc_spool_piece_correctors(lhc, twiss_b1=tw_b1, twiss_b2=None)

# k1s local + global correction (uses MQS)
xmlhc.correct_k1s(lhc, twiss_b1=tw_b1, twiss_b2=None)

# k2s local + global correction (uses MSS)
xmlhc.correct_k2s(lhc, twiss_b1=tw_b1, twiss_b2=None)

line.to_json('hllhc_b1_v19_round_imo300_with_arc_errors_and_corrections.json')