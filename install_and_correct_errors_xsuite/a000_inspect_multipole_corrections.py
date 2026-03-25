import xtrack as xt
import xmask as xm
import numpy as np

line_test = xt.load("hllhc_b1_v19_round_imo300_with_arc_errors_and_corrections.json")

env_test = line_test.env
env_test['b1'] = line_test
line_name = 'b1'

# Chromatic coupling integrand
tw_test = line_test.twiss4d(strengths=True)

tw_test['chrom_coupl_source'] = (tw_test.k2sl * tw_test.dx * np.sqrt(tw_test.betx * tw_test.bety)
                            * np.exp(1j*2*np.pi*(tw_test.mux - tw_test.muy)))

# Chromatic coupling integral arc by arc
chrom_coupling_integ_test = {}
chrom_coupling_integ_ref = {}
for arc in ['12', '23', '34', '45', '56', '67', '78', '81']:
    if line_name == 'b1':
        start = f's.ds.r{arc[0]}.b1'
        end = f'e.ds.l{arc[1]}.b1'
    else:
        start = f'e.ds.l{arc[1]}.b2'
        end = f's.ds.r{arc[0]}.b2'
    tt_test_arc = tw_test.rows[start:end]

    chrom_coupling_integ_test[arc] = np.sum(tt_test_arc.chrom_coupl_source)

twom_test = line_test.twiss4d(coupling_edw_teng=True, delta0=0.5e-3)

nlchr_test = line_test.get_non_linear_chromaticity(delta0_range=(-5e-4, 5e-4), num_delta=20)

cminus_delta_test = np.array([ttt.c_minus for ttt in nlchr_test.twiss])

# Check on local chromatic coupling integrals
arc_names = list(chrom_coupling_integ_test.keys())
arc_chrom_coupl_test = [chrom_coupling_integ_test[arc_name] for arc_name in arc_names]

global_chrom_coupling = {}
global_chrom_coupling_no_corr = {}
local_chrom_coupling = {}

# Measure chromatic coupling without correction
tt_on_corr = env_test.vars.get_table().rows['on_corr_k2s.*']
env_test.set(tt_on_corr.name, 0)
nlchr_test_no_corr = line_test.get_non_linear_chromaticity(num_delta=20, delta0_range=(-5e-4, 5e-4))
env_test.set(tt_on_corr.name, 1)
cminus_delta_test_no_corr = np.array([ttt.c_minus for ttt in nlchr_test_no_corr.twiss])

global_chrom_coupling[line_name] = cminus_delta_test
local_chrom_coupling[line_name] = arc_chrom_coupl_test
global_chrom_coupling_no_corr[line_name] = cminus_delta_test_no_corr

assert np.max(global_chrom_coupling_no_corr['b1']) > 2e-3
assert np.all(np.abs(global_chrom_coupling['b1']) < 4e-4)

# Check effect of mcs on chromaticity
env_test.set(env_test.vars.get_table().rows['on_corr_kcs.*'], 0)
tw_b1_mcs_off = env_test['b1'].twiss4d()
env_test.set(env_test.vars.get_table().rows['on_corr_kcs.*'], 1)
tw_b1_mcs_on = env_test['b1'].twiss4d()

assert np.abs(tw_b1_mcs_off.dqx) > 100
assert np.abs(tw_b1_mcs_off.dqy) > 100
assert np.abs(tw_b1_mcs_on.dqx) < 20
assert np.abs(tw_b1_mcs_on.dqy) < 20

import matplotlib.pyplot as plt
plt.close('all')

line_name = 'b1'

c_minus_delta_test = global_chrom_coupling[line_name]
c_minus_delta_test_no_corr = global_chrom_coupling_no_corr[line_name]
arc_chrom_coupl_test = local_chrom_coupling[line_name]

plt.figure(1)
plt.plot(1e3 * nlchr_test_no_corr.delta0, c_minus_delta_test_no_corr, label='Test no corr', linestyle='dashed')
plt.plot(1e3 * nlchr_test.delta0, c_minus_delta_test, label='Test')
plt.xlabel(r'$\delta_0$ [10$^{-3}$]')
plt.ylabel('C-')
plt.legend()
plt.suptitle(line_name)

# Bar plot of chromatic coupling integral arc by arc
plt.figure(2)

x = np.arange(len(arc_names))
plt.bar(x - 0.2, np.abs(arc_chrom_coupl_test), width=0.4, label='Test')
plt.xticks(x, arc_names)
plt.xlabel('Arc name')
plt.ylabel('Abs of chromatic coupling integral')
plt.legend()
plt.suptitle(line_name)
