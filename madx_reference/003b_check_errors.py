import xtrack as xt
import xobjects as xo
import numpy as np

line_test = xt.load('test_line_with_errors.json')
line_ref = xt.load('hllhc_b1_v19_round_imo300_arcs.json')

tt_test = line_test.get_table(attr=True)
tt_ref = line_ref.get_table(attr=True)

for arc in ['12', '23', '34', '45', '56', '67', '78', '81']:
    start = f's.ds.r{arc[0]}.b1'
    end = f's.ds.r{arc[1]}.b1'
    tt_test_arc = tt_test.rows[f'{start}:{end}']
    tt_ref_arc = tt_ref.rows[f'{start}:{end}']
    for ii in range(2, 15): # min_order to max_order
        xo.assert_allclose(tt_test_arc[f'k{ii}l'], tt_ref_arc[f'k{ii}l'], rtol=1e-10, atol=1e-10)
        xo.assert_allclose(tt_test_arc[f'k{ii}sl'], tt_ref_arc[f'k{ii}sl'], rtol=1e-10, atol=1e-10)
