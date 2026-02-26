import numpy as np
import re

import xtrack as xt

def assign_spool_pieces(line, beam):
    arc_map = {
        ('r', '1'): 12, ('l', '2'): 12,
        ('r', '2'): 23, ('l', '3'): 23,
        ('r', '3'): 34, ('l', '4'): 34,
        ('r', '4'): 45, ('l', '5'): 45,
        ('r', '5'): 56, ('l', '6'): 56,
        ('r', '6'): 67, ('l', '7'): 67,
        ('r', '7'): 78, ('l', '8'): 78,
        ('r', '8'): 81, ('l', '1'): 81
    }
    mb_pattern = re.compile(r"mb\.\w*([rl])([1-8])\.[\w\.]+")
    mcs_pattern = re.compile(r"mcs\.\w*([rl])([1-8])\.[\w\.]+")
    mco_pattern = re.compile(r"mco\.\w*([rl])([1-8])\.[\w\.]+")
    mcd_pattern = re.compile(r"mcd\.\w*([rl])([1-8])\.[\w\.]+")

    mb_per_arc = {
        12: [], 
        23: [], 
        34: [], 
        45: [], 
        56: [], 
        67: [], 
        78: [], 
        81: [], 
    }
    mcs_per_arc = {
        12: [], 
        23: [], 
        34: [], 
        45: [], 
        56: [], 
        67: [], 
        78: [], 
        81: [], 
    }
    mco_per_arc = {
        12: [], 
        23: [], 
        34: [], 
        45: [], 
        56: [], 
        67: [], 
        78: [], 
        81: [], 
    }
    mcd_per_arc = {
        12: [], 
        23: [], 
        34: [], 
        45: [], 
        56: [], 
        67: [], 
        78: [], 
        81: [], 
    }

    for name in line.element_names:
        match_mb = mb_pattern.fullmatch(name)
        match_mcs = mcs_pattern.fullmatch(name)
        match_mco = mco_pattern.fullmatch(name)
        match_mcd = mcd_pattern.fullmatch(name)
        
        if match_mb:
            letter, digit = match_mb.groups()
            arc = arc_map.get((letter, digit))
            if arc:
                ele_type = line[name].to_dict()["__class__"]
                if ele_type == "Drift" or ele_type == "Marker" or ele_type.startswith("Limit"):
                    continue
                mb_per_arc[arc].append(name)
        elif match_mcs:
            letter, digit = match_mcs.groups()
            arc = arc_map.get((letter, digit))
            if arc:
                ele_type = line[name].to_dict()["__class__"]
                if ele_type == "Drift" or ele_type == "Marker" or ele_type.startswith("Limit"):
                    continue
                mcs_per_arc[arc].append(name)
        elif match_mco:
            letter, digit = match_mco.groups()
            arc = arc_map.get((letter, digit))
            if arc:
                ele_type = line[name].to_dict()["__class__"]
                if ele_type == "Drift" or ele_type == "Marker" or ele_type.startswith("Limit"):
                    continue
                mco_per_arc[arc].append(name)
        elif match_mco:
            letter, digit = match_mco.groups()
            arc = arc_map.get((letter, digit))
            if arc:
                ele_type = line[name].to_dict()["__class__"]
                if ele_type == "Drift" or ele_type == "Marker" or ele_type.startswith("Limit"):
                    continue
                mco_per_arc[arc].append(name)
        elif match_mcd:
            letter, digit = match_mcd.groups()
            arc = arc_map.get((letter, digit))
            if arc:
                ele_type = line[name].to_dict()["__class__"]
                if ele_type == "Drift" or ele_type == "Marker" or ele_type.startswith("Limit"):
                    continue
                mcd_per_arc[arc].append(name)
        else:
            continue

    for arc in [12, 23, 34, 45, 56, 67, 78, 81]:
        total_k2l = 0.0
        total_k3l = 0.0
        total_k4l = 0.0

        for name in mb_per_arc[arc]:
            total_k2l += line[name].knl[2]
            total_k3l += line[name].knl[3]
            total_k4l += line[name].knl[4]

        k2l_corr = total_k2l / len(mcs_per_arc[arc]) if len(mcs_per_arc[arc]) != 0 else 0 
        k3l_corr = total_k3l / len(mco_per_arc[arc]) if len(mco_per_arc[arc]) != 0 else 0
        k4l_corr = total_k4l / len(mcd_per_arc[arc]) if len(mcd_per_arc[arc]) != 0 else 0

        if beam == 1:
            kcs = -k2l_corr / line.vv["l.mcs"]
            kco = -k3l_corr / line.vv["l.mco"]
            kcd = -k4l_corr / line.vv["l.mcd"]
        else:
            kcs = k2l_corr / line.vv["l.mcs"]
            kco = k3l_corr / line.vv["l.mco"]
            kcd = k4l_corr / line.vv["l.mcd"]

        if np.isfinite(kcs):
            line.vv[f"kcs.a{arc}b{beam}"] = kcs
        if np.isfinite(kco):
            line.vv[f"kco.a{arc}b{beam}"] = kco
        if np.isfinite(kcd):
            line.vv[f"kcd.a{arc}b{beam}"] = kcd


def save_crossing(line):
    crossing_knobs = ['on_a2', 'on_a8', 'on_disp', 'on_o2', 'on_o8', 'on_oh1', 'on_oh2', 'on_oh5', 'on_oh8', 'on_ov1',
                    'on_ov2', 'on_ov5', 'on_ov8', 'on_sep1', 'on_sep2h', 'on_sep2v', 'on_sep5', 'on_sep8h', 'on_sep8v',
                    'on_ssep1', 'on_ssep5', 'on_x1', 'on_x2h', 'on_x2v', 'on_x5', 'on_x8h', 'on_x8v', 'on_xx1', 'on_xx5', 
                    'on_alice', 'on_lhcb']
    for var in crossing_knobs:
        line.vars[var + "_stored"] = line.vv[var]


def disable_crossing(line):
    crossing_knobs = ['on_a2', 'on_a8', 'on_disp', 'on_o2', 'on_o8', 'on_oh1', 'on_oh2', 'on_oh5', 'on_oh8', 'on_ov1',
                    'on_ov2', 'on_ov5', 'on_ov8', 'on_sep1', 'on_sep2h', 'on_sep2v', 'on_sep5', 'on_sep8h', 'on_sep8v',
                    'on_ssep1', 'on_ssep5', 'on_x1', 'on_x2h', 'on_x2v', 'on_x5', 'on_x8h', 'on_x8v', 'on_xx1', 'on_xx5', 
                    'on_alice', 'on_lhcb']
    for var in crossing_knobs:
        line.vv[var] = 0.0


def enable_crossing(line):
    crossing_knobs = ['on_a2', 'on_a8', 'on_disp', 'on_o2', 'on_o8', 'on_oh1', 'on_oh2', 'on_oh5', 'on_oh8', 'on_ov1',
                    'on_ov2', 'on_ov5', 'on_ov8', 'on_sep1', 'on_sep2h', 'on_sep2v', 'on_sep5', 'on_sep8h', 'on_sep8v',
                    'on_ssep1', 'on_ssep5', 'on_x1', 'on_x2h', 'on_x2v', 'on_x5', 'on_x8h', 'on_x8v', 'on_xx1', 'on_xx5', 
                    'on_alice', 'on_lhcb']
    for var in crossing_knobs:
        line.vv[var] = line.vv[var + "_stored"]


def set_correctors(line):
    crossing_knobs = {'on_a2', 'on_a8', 'on_disp', 'on_o2', 'on_o8', 'on_oh1', 'on_oh2', 'on_oh5', 'on_oh8', 'on_ov1',
                    'on_ov2', 'on_ov5', 'on_ov8', 'on_sep1', 'on_sep2h', 'on_sep2v', 'on_sep5', 'on_sep8h', 'on_sep8v',
                    'on_ssep1', 'on_ssep5', 'on_x1', 'on_x2h', 'on_x2v', 'on_x5', 'on_x8h', 'on_x8v', 'on_xx1', 'on_xx5'}
    crossing_currents = set()
    for nn in line.vars.get_table().name:
        if line.vars[nn]._expr is None:
            continue
        deps = {vv._key for vv in line.vars[nn]._expr._get_dependencies()}
        if any([vvv in crossing_knobs for vvv in deps]):
            crossing_currents.add(nn)
    crossing_correctors = {vvv._key for vv in crossing_currents for vvv in line.vars[vv]._find_dependant_targets()}

    tt = line.get_table()
    mask = ~np.array([nn.startswith('Limit') or nn.startswith('Drift') or nn.startswith('Marker')
                        for nn in tt.element_type])
    tt_h_correctors = tt.rows[mask].rows['mcb.*'].rows['.*h\..*'].name
    line.steering_correctors_x = list({nn for nn in tt_h_correctors if nn not in crossing_correctors})
    tt_v_correctors = tt.rows[mask].rows['mcb.*'].rows['.*v\..*'].name
    line.steering_correctors_y = list({nn for nn in tt_v_correctors if nn not in crossing_correctors})

    mask = ~np.array([nn.startswith('Limit') for nn in tt.element_type])
    tt_monitors = tt.rows[mask].rows['bpm\..*'].rows['.*(?<!_entry)$'].rows['.*(?<!_exit)$'].name
    tt_monitors = [nn for nn in tt_monitors if not nn.startswith('bpmwa')]
    tt_monitors = [nn for nn in tt_monitors if not nn.startswith('bpmwb')]
    tt_monitors = [nn for nn in tt_monitors if not nn.startswith('bpmse')]
    tt_monitors = [nn for nn in tt_monitors if not nn.startswith('bpmsd')]
    # Remove double-counted monitors (multiple at same s), keep shortest name
    monitors_by_s = {}
    for nn, ss in zip(tt.name, tt.s):
        if nn in tt_monitors:
            ss_key = [sss for sss in monitors_by_s.keys() if np.isclose(ss, sss, atol=0.001)]
            if not ss_key:
                monitors_by_s[ss] = [nn]
            else:
                monitors_by_s[ss_key[0]].append(nn)
    tt_monitors = [min(nn, key=len) for nn in monitors_by_s.values()]
    line.steering_monitors_x = tt_monitors
    line.steering_monitors_y = tt_monitors


def match_tune_chrom(line, qx, qy, dqx, dqy, tol=1e-3):
    if hasattr(tol, '__iter__'):
        penalty = 1e10
        for this_tol in tol:
            if penalty < this_tol:
                continue
            opt = match_tune_chrom(line, qx=qx, qy=qy, dqx=dqx, dqy=dqy,
                                   tol=this_tol)
            if this_tol != tol[-1]:
                penalty = np.sqrt((opt.target_status(True).residue**2).sum())
        return opt
    return line.match(
        method='6d', # <- passed to twiss
        vary=[
            xt.VaryList(['kqtf', 'kqtd'], step=tol*1e-1, tag='quad'),
            xt.VaryList(['ksf', 'ksd'], step=tol*1e-1, tag='sext'),
        ],
        targets = [
            xt.TargetSet(qx=qx, qy=qy, tol=tol, tag='tune'),
            xt.TargetSet(dqx=dqx, dqy=dqy, tol=tol, tag='chrom'),
        ]
    )


def match_coupling(line, beam, c_minus, tol=1e-3):
    if hasattr(tol, '__iter__'):
        penalty = 1e10
        for this_tol in tol:
            if penalty < this_tol:
                continue
            opt = match_coupling(line, c_minus=c_minus, tol=this_tol)
            if this_tol != tol[-1]:
                penalty = np.sqrt((opt.target_status(True).residue**2).sum())
        return opt
    return line.match(
        method='6d',
        vary=[xt.VaryList([f'cmrskew', f'cmiskew'], limits=[-0.5e-2, 0.5e-2], step=tol*1e-1)],
        targets=[
            xt.Target('c_minus_re_0', c_minus, tol=tol), xt.Target('c_minus_im_0', 0, tol=tol)]
    )