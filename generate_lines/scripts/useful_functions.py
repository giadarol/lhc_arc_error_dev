"""
Unified helper functions for HL-LHC error line generation.

Consolidates:
    - Mask running / symlinks / cleanup (from errors_in_arcs/functions.py)
    - D1 table generation (from errors_in_D1/d1_errors.py + d1_table_functions.py)
    - D2 table generation (from errors_in_D2/d2_errors.py + d2_table_functions.py)
    - Triplet table generation (from errors_in_triplet/triplet_errors.py + triplet_table_functions.py)
    - Mask content generation (new — builds MAD-X mask string from config)
"""

import json
import os
import shutil
from pathlib import Path

import numpy as np
import yaml
from cpymad.madx import Madx
import xtrack as xt


# =============================================================================
# Constants
# =============================================================================

TRACKING_TOOLS = "./scripts/madx"
TRACKING_DATA = "./data"
MAX_ORDER = 15  # lhcerrors tables always have orders 1..15
POSITIONS = ["r5", "l5", "r1", "l1"]  # standard iteration order for triplet


# =============================================================================
# Symlink management
# =============================================================================

def setup_symlinks(working_dir: Path = Path("."), acc_models_lhc: str = None):
    """Create symlinks required by MAD-X submodules (errors, tools, modules, acc-models-lhc)."""
    links = {
        "errors":  f"{TRACKING_TOOLS}/errors",
        "tools":   f"{TRACKING_TOOLS}/tools",
        "modules": f"{TRACKING_TOOLS}/modules",
        "madxscripts": f"{TRACKING_TOOLS}",
        "madxdata": f"{TRACKING_DATA}/madx",
    }
    if acc_models_lhc is not None:
        links["acc-models-lhc"] = acc_models_lhc
    for name, target in links.items():
        link_path = working_dir / name
        if link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(target)


def cleanup(working_dir: Path = Path(".")):
    """Remove temporary files, symlinks and dirs created during a mask run."""
    temp_files = [
        "params_generated.madx",
        "error_all.tfs",
        "twiss_lhcb1.tfs",
        "twiss_lhcb2.tfs",
        "twiss_lhcb4.tfs",
        "twiss_crossing_*",
        "internal_mag_pot.tfs",
        "mad_*.tfs",
        "fc.*",
        "last_twiss.*.gz",
        "d2_aperture",
    ]
    temp_symlinks = ["errors", "tools", "modules", "acc-models-lhc", "madxdata", "madxscripts"]
    temp_dirs = ["temp"]

    for pattern in temp_files:
        if "*" in pattern:
            for f in working_dir.glob(pattern):
                f.unlink()
        else:
            f = working_dir / pattern
            if f.exists() and not f.is_symlink():
                f.unlink()

    for name in temp_symlinks:
        link = working_dir / name
        if link.is_symlink():
            link.unlink()

    for d in temp_dirs:
        dp = working_dir / d
        if dp.exists():
            shutil.rmtree(dp)


# =============================================================================
# Mask execution
# =============================================================================

def run_madx_mask(mask_file, params=None, working_dir: Path = Path("."), verbose: bool = True,
                  acc_models_lhc: str = None):
    """Run a MAD-X mask file via cpymad.

    Creates symlinks, sets MAD-X variables from *params*, then calls the mask.
    Returns the Madx instance.
    """
    setup_symlinks(working_dir, acc_models_lhc=acc_models_lhc)

    original_dir = os.getcwd()
    os.chdir(working_dir)
    try:
        madx = Madx(stdout=verbose)
        if params:
            for name, value in params.items():
                madx.input(f"{name} = {value};")
        if not os.path.exists(str(mask_file)):
            print("Mask file does not exist!")
        madx.call(str(mask_file))
        return madx
    except Exception as e:
        print(f"\nERROR during MAD-X execution: {e}")
        raise
    finally:
        os.chdir(original_dir)


# =============================================================================
# xtrack conversion
# =============================================================================

def madx_to_xtrack_line(madx, sequence_name: str):
    """Convert a MAD-X sequence to an xtrack Line."""
    beam = madx.sequence[sequence_name].beam

    line = xt.Line.from_madx_sequence(
        madx.sequence[sequence_name],
        deferred_expressions=True,
        install_apertures=True,
        apply_madx_errors=True,
    )
    line.particle_ref = xt.Particles(
        mass0=beam.mass * 1e9,
        q0=beam.charge,
        energy0=beam.energy * 1e9,
    )
    return line


# =============================================================================
# D1 (MBXF / MBXAB) — FQ loading and table generation
# =============================================================================

def load_d1_fq(magnet_name, fq_dir, measurement="cold_nominal"):
    """Load FQ data for a single MBXF magnet. Returns {n: (bn, an)}."""
    fq_dir = Path(fq_dir)
    path = fq_dir / f"FQ_{magnet_name}_{measurement}.json"
    if not path.exists():
        raise FileNotFoundError(f"No FQ file for {magnet_name}: {path}")
    with open(path) as f:
        data = json.load(f)
    fq = {}
    for mult in data["multipoles"]:
        n = mult["n"]
        fq[n] = (mult.get("bn", 0.0), mult.get("an", 0.0))
    return fq


def get_d1_is_inv(position):
    """D1 is_inv: L-side -> 0, R-side -> 1."""
    side = position[1].upper()
    return 0 if side == "L" else 1


def _write_d1_error_table(filepath, fq_data, magnet_name, label):
    """Write one D1 MAD-X error table file (M=measured, U=R=0)."""
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    s = "_MBXAB_"
    lines = [f"! Measured FQ for MBXAB.{label} (magnet: {magnet_name})"]
    for ba, idx in [("b", 0), ("a", 1)]:
        for n in range(1, MAX_ORDER + 1):
            val = fq_data[n][idx] if n in fq_data else 0.0
            lines.append(
                f"{ba}{n}M{s}col := {val: .15g} ; "
                f"{ba}{n}U{s}col := 0.0 ; "
                f"{ba}{n}R{s}col := 0.0 ;"
            )
    for ba in ("b", "a"):
        for n in range(1, MAX_ORDER + 1):
            lines.append(
                f"{ba}{n}M{s}inj := {ba}{n}M{s}col ; "
                f"{ba}{n}U{s}inj := 0.0 ; "
                f"{ba}{n}R{s}inj := 0.0 ;"
            )
    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")


def build_d1_tables(config, output_dir="tables/D1_tables"):
    """Generate all D1 error table files."""
    output_dir = Path(output_dir)
    d1 = config["D1"]
    fq_dir = d1["fq_path"]
    measurement = d1["measurement"]
    positions = d1["positions"]

    for pos, magnet_name in positions.items():
        label = pos.upper()
        table_file = f"{output_dir}/MBXAB_{label}_errortable"

        fq = load_d1_fq(magnet_name, fq_dir, measurement)
        _write_d1_error_table(table_file, fq, magnet_name, label)
        print(f"  MBXAB.{label}: {magnet_name}")


# =============================================================================
# D2 (MBRD) — FQ loading and table generation
# =============================================================================

def load_d2_fq(magnet_name, aperture, fq_dir, measurement="cold_nominal_extrapolation"):
    """Load FQ data for a single MBRD magnet aperture. Returns {n: (bn, an)}."""
    fq_dir = Path(fq_dir)
    path = fq_dir / f"FQ_{magnet_name}_{aperture}_{measurement}.json"
    if not path.exists():
        raise FileNotFoundError(f"No FQ file for {magnet_name} {aperture}: {path}")
    with open(path) as f:
        data = json.load(f)
    fq = {}
    for mult in data["multipoles"]:
        n = mult["n"]
        fq[n] = (mult.get("bn", 0.0), mult.get("an", 0.0))
    return fq


def get_d2_is_inv(position):
    """D2 is_inv: L-side -> 1, R-side -> 0 (opposite to D1)."""
    side = position[1].upper()
    return 1 if side == "L" else 0


def _write_d2_error_table(filepath, fq_ap1, fq_ap2, magnet_name, label):
    """Write one D2 MAD-X error table file with V1/V2 columns (M=measured, U=R=0)."""
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    s = "_MBRD_"
    lines = [f"! Measured FQ for MBRD.{label} (magnet: {magnet_name}, AP1->V1, AP2->V2)"]
    for ba, idx in [("b", 0), ("a", 1)]:
        for n in range(1, MAX_ORDER + 1):
            v1 = fq_ap1[n][idx] if n in fq_ap1 else 0.0
            v2 = fq_ap2[n][idx] if n in fq_ap2 else 0.0
            lines.append(
                f"{ba}{n}M{s}V1_col := {v1: .15g} ; "
                f"{ba}{n}M{s}V2_col := {v2: .15g} ; "
                f"{ba}{n}U{s}col := 0.0 ; "
                f"{ba}{n}R{s}col := 0.0 ;"
            )
    for ba in ("b", "a"):
        for n in range(1, MAX_ORDER + 1):
            lines.append(
                f"{ba}{n}M{s}V1_inj := {ba}{n}M{s}V1_col ; "
                f"{ba}{n}M{s}V2_inj := {ba}{n}M{s}V2_col ; "
                f"{ba}{n}U{s}inj := 0.0 ; "
                f"{ba}{n}R{s}inj := 0.0 ;"
            )
    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")


def build_d2_tables(config, output_dir="tables/D2_tables"):
    """Generate all D2 error table files."""
    output_dir = Path(output_dir)
    d2 = config["D2"]
    fq_dir = d2["fq_path"]
    measurement = d2["measurement"]
    positions = d2["positions"]

    beam1_ap = d2.get("beam1_aperture", {})

    for pos, magnet_name in positions.items():
        label = pos.upper()
        table_file = f"{output_dir}/MBRD_{label}_errortable"

        fq_ap1 = load_d2_fq(magnet_name, "AP1", fq_dir, measurement)
        fq_ap2 = load_d2_fq(magnet_name, "AP2", fq_dir, measurement)

        side = "left" if pos[1].upper() == "L" else "right"
        b1_ap = beam1_ap.get(side, "AP1")
        if b1_ap == "AP2":
            fq_v1, fq_v2 = fq_ap2, fq_ap1
        else:
            fq_v1, fq_v2 = fq_ap1, fq_ap2

        _write_d2_error_table(table_file, fq_v1, fq_v2, magnet_name, label)
        print(f"  MBRD.{label}: {magnet_name} ({b1_ap}->V1, {'AP2' if b1_ap == 'AP1' else 'AP1'}->V2)")


# =============================================================================
# Triplet (MQXFA / MQXFB) — FQ loading, synthetic generation, table generation
# =============================================================================

def load_triplet_fq(magnet_name, fq_dir, measurement="cold_nominal"):
    """Load FQ data for a single MQXFA/MQXFB magnet. Returns {n: (bn, an)}.

    Handles both naming conventions:
      - MQXFA: {magnet}_{CA}.json (search with glob)
      - MQXFB: FQ_{magnet}_{measurement}.json
    """
    fq_dir = Path(fq_dir)

    # Try MQXFB naming first
    mqxfb_path = fq_dir / f"FQ_{magnet_name}_{measurement}.json"
    if mqxfb_path.exists():
        path = mqxfb_path
    else:
        # MQXFA naming: {magnet}_*.json (CA suffix varies)
        matches = list(fq_dir.glob(f"{magnet_name}_*.json"))
        if len(matches) == 0:
            raise FileNotFoundError(
                f"No FQ file found for {magnet_name} in {fq_dir}"
            )
        if len(matches) > 1:
            raise FileNotFoundError(
                f"Multiple FQ files for {magnet_name} in {fq_dir}: {matches}"
            )
        path = matches[0]

    with open(path) as f:
        data = json.load(f)

    fq = {}
    for mult in data["multipoles"]:
        n = mult["n"]
        fq[n] = (mult.get("bn", 0.0), mult.get("an", 0.0))
    return fq


def _sample_truncated_normal(mean, sigma, trunc_sigma, size, rng):
    """Sample from N(mean, sigma) truncated at +/-trunc_sigma standard deviations."""
    samples = np.empty(size)
    count = 0
    while count < size:
        batch = rng.normal(mean, sigma, size=size - count)
        mask = np.abs(batch - mean) <= trunc_sigma * sigma
        n_good = mask.sum()
        samples[count:count + n_good] = batch[mask]
        count += n_good
    return samples


def generate_synthetic_magnets(random_model_path, magnet_names, seed=42):
    """Generate synthetic FQ data from the LHC statistical model.

    Model: b_n = S + (xi_U/1.5)*U + xi_R*R
      - xi_U ~ N(0,1) truncated +/-1.5 (common per multipole across all magnets)
      - xi_R ~ N(0,1) truncated +/-3   (independent per magnet)
      - R = (T - U) / 3

    Returns dict {magnet_name: {n: (bn, an)}}.
    """
    with open(random_model_path) as f:
        model = yaml.safe_load(f)

    mm = model["multipole_model"]
    orders = mm["orders"]
    components = mm["components"]

    rng = np.random.default_rng(seed)

    # Draw xi_U (one per multipole order, shared across magnets)
    n_components = len(orders) * 2  # bn + an for each order
    xi_U = _sample_truncated_normal(0, 1, 1.5, n_components, rng)

    result = {}
    for i_mag, mag_name in enumerate(magnet_names):
        fq = {}
        for i_ord, n in enumerate(orders):
            bn_params = components[f"b{n}"]
            an_params = components[f"a{n}"]

            S_b, U_b, T_b = bn_params["S"], bn_params["U"], bn_params["T"]
            R_b = (T_b - U_b) / 3.0
            xi_R_b = _sample_truncated_normal(0, 1, 3.0, 1, rng)[0]
            bn_val = S_b + (xi_U[2 * i_ord] / 1.5) * U_b + xi_R_b * R_b

            S_a, U_a, T_a = an_params["S"], an_params["U"], an_params["T"]
            R_a = (T_a - U_a) / 3.0
            xi_R_a = _sample_truncated_normal(0, 1, 3.0, 1, rng)[0]
            an_val = S_a + (xi_U[2 * i_ord + 1] / 1.5) * U_a + xi_R_a * R_a

            fq[n] = (bn_val, an_val)

        result[mag_name] = fq

    return result


def get_triplet_is_inv(half, position, q_type):
    """Determine is_inv for a triplet half/position.

    Q1/Q3: is_inv = (half=='a') != (side=='R')
    Q2:    is_inv = (half=='a') == (side=='R')
    """
    side_is_right = position[0].upper() == "R"
    half_is_a = (half == "a")

    if q_type in ("Q1", "Q3"):
        is_inv = int(half_is_a != side_is_right)
    else:
        is_inv = int(half_is_a == side_is_right)

    return is_inv


def _write_triplet_error_table(filepath, magnet_type, fq_data, magnet_name, label):
    """Write one triplet MAD-X error table file (M=measured, U=R=0)."""
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    s = f"_{magnet_type}_"
    lines = [f"! Measured FQ for {magnet_type}.{label} (magnet: {magnet_name})"]
    for ba, idx in [("b", 0), ("a", 1)]:
        for n in range(1, MAX_ORDER + 1):
            val = fq_data[n][idx] if n in fq_data else 0.0
            lines.append(
                f"{ba}{n}M{s}col := {val: .15g} ; "
                f"{ba}{n}U{s}col := 0.0 ; "
                f"{ba}{n}R{s}col := 0.0 ;"
            )
    for ba in ("b", "a"):
        for n in range(1, MAX_ORDER + 1):
            lines.append(
                f"{ba}{n}M{s}inj := {ba}{n}M{s}col ; "
                f"{ba}{n}U{s}inj := 0.0 ; "
                f"{ba}{n}R{s}inj := 0.0 ;"
            )
    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")

def build_triplet_tables(config, output_dir="tables/IT_tables", synthetic_seed=None):
    """Generate all IT error table files."""
    output_dir = Path(output_dir)
    triplet = config["triplet"]
    cas = triplet["cryoassemblies"]
    fq_a = triplet["fq_paths"]["MQXFA"]
    fq_b = triplet["fq_paths"]["MQXFB"]
    measurement = triplet["measurement"]

    synth_config = triplet["synthetic"]
    seed = synthetic_seed if synthetic_seed is not None else synth_config["seed"]
    synth_fq = generate_synthetic_magnets(
        triplet["fq_paths"]["random_model"], synth_config["magnets"], seed=seed,
    )
    synth_set = set(synth_config["magnets"])

    def load_magnet_fq(name):
        if name in synth_set:
            return synth_fq[name]
        if name.upper().startswith("MQXFB"):
            return load_triplet_fq(name, fq_b, measurement)
        return load_triplet_fq(name, fq_a, measurement)

    for q_type in ("Q1", "Q2", "Q3"):
        q_config = triplet[q_type]
        mt = q_config["prefix"].upper()
        slice_ids = q_config["slice_ids"]

        for pos in POSITIONS:
            assignment = q_config["positions"][pos]
            magnets = cas[assignment] if isinstance(assignment, str) else assignment

            for mag_name, sid in zip(magnets, slice_ids):
                label = f"{sid}{pos}".upper()
                table_file = f"{output_dir}/{mt}_{label}_errortable"

                fq = load_magnet_fq(mag_name)
                _write_triplet_error_table(table_file, mt, fq, mag_name, label)
                print(f"  {mt}.{label}: {mag_name}")

def _write_mcbrd_error_table(filepath, mcbrdh_v1, mcbrdh_v2, mcbrdv_v1, mcbrdv_v2, label):
    """Write one per-position MCBRD MAD-X error table file with MCBRDH + MCBRDV sections."""
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    lines = [f"! MCBRD error table for position {label}"]

    for mtype, v1_data, v2_data in [("MCBRDH", mcbrdh_v1, mcbrdh_v2),
                                     ("MCBRDV", mcbrdv_v1, mcbrdv_v2)]:
        s = f"_{mtype}_"
        for ba, idx in [("b", 0), ("a", 1)]:
            for n in range(1, MAX_ORDER + 1):
                v1 = v1_data.get(n, (0.0, 0.0))[idx]
                v2 = v2_data.get(n, (0.0, 0.0))[idx]
                lines.append(
                    f"{ba}{n}M{s}V1_col := {v1: .6f} ; "
                    f"{ba}{n}M{s}V2_col := {v2: .6f} ; "
                    f"{ba}{n}U{s}col := 0.0 ; "
                    f"{ba}{n}R{s}col := 0.0 ;"
                )
        for ba, idx in [("b", 0), ("a", 1)]:
            for n in range(1, MAX_ORDER + 1):
                lines.append(
                    f"{ba}{n}M{s}V1_inj := {ba}{n}M{s}V1_col ; "
                    f"{ba}{n}M{s}V2_inj := {ba}{n}M{s}V2_col ; "
                    f"{ba}{n}U{s}inj := 0.0 ; "
                    f"{ba}{n}R{s}inj := 0.0 ;"
                )

    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")

def build_mcbrd_tables(config, output_dir="tables/MCBRD_tables", seed=1):
    """Generate per-position MCBRD error tables.

    Uniform random multipole errors centered at 0 (no systematic component).
    BRef = 2 Tm. Error ranges: +-3 units for all orders, except b2/a2/b3/a3: +-15 units.
    """
    output_dir = Path(output_dir)
    mcbrd = config["MCBRD"]
    positions = mcbrd["positions"]

    default_range = mcbrd.get("default_range", 3.0)
    special_range = mcbrd.get("special_range", 15.0)
    special_orders = mcbrd.get("special_orders", {2, 3})

    rng = np.random.default_rng(seed)

    def random_multipoles_uniform():
        v1, v2 = {}, {}
        for n in range(1, MAX_ORDER + 1):
            half_range = special_range if n in special_orders else default_range
            v1[n] = (rng.uniform(-half_range, half_range),
                     rng.uniform(-half_range, half_range))
            v2[n] = (rng.uniform(-half_range, half_range),
                     rng.uniform(-half_range, half_range))
        return v1, v2

    for pos in positions:
        h_v1, h_v2 = random_multipoles_uniform()
        v_v1, v_v2 = random_multipoles_uniform()
        table_file = output_dir / f"MCBRD_{pos}_errortable"
        _write_mcbrd_error_table(table_file, h_v1, h_v2, v_v1, v_v2, pos)
        print(f"  MCBRD.{pos}: table written")


