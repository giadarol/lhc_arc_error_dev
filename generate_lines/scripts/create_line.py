import os
import yaml
import argparse

import xobjects as xo
import xdeps as xd
import xtrack as xt
import xpart as xp
import xfields as xf
import xcoll as xc
import xsuite as xs

from correction_tools import *
from useful_functions import (
    build_d1_tables,
    build_d2_tables,
    build_mcbrd_tables,
    build_triplet_tables,
    run_madx_mask,
    madx_to_xtrack_line,
    cleanup,
)

"""
Generates an xtrack line with any combination of error sources:
  - Arc errors (random, seed-dependent)
  - D1 measured FQ (MBXF)
  - D2 measured FQ (MBRD, 2-in-1)
  - Triplet measured FQ (MQXFA/MQXFB)

Run cell-by-cell in an interactive environment (VS Code / Jupyter).
Launch from the making_error_line/HL_FQ_lattice/ directory.
"""

parser = argparse.ArgumentParser()
parser.add_argument(
    "beam", 
    type=int, 
    help="1 or 2"
)
parser.add_argument(
    "seed", 
    type=int, 
    help="Error seed."
)
parser.add_argument(
    "version", 
    type=int, 
    help="HL-LHC lattice version, e.g. '16' or '19'."
)
parser.add_argument(
    "optics", 
    type=str, 
    help="HL-LHC optics, can be 'round' or 'flat'."
)
parser.add_argument(
    "i_mo", 
    type=float, 
    help="Octupole current in A."
)
parser.add_argument(
    "ref_line_path", 
    type=str, 
    help="Path to reference line without errors."
)
parser.add_argument(
    "out_path", 
    type=str, 
    help="Output path."
)


def main():
    '''
    Print xsuite package versions
    '''
    print(f"xobjects: {xo.__version__}")
    print(f"xdeps: {xd.__version__}")
    print(f"xtrack: {xt.__version__}")
    print(f"xpart: {xp.__version__}")
    print(f"xfields: {xf.__version__}")
    print(f"xcoll: {xc.__version__}")
    print(f"xsuite: {xs.__version__}")


    '''
    Parse arguments
    '''
    args = parser.parse_args()
    beam = args.beam
    seed = args.seed
    version = args.version
    optics = args.optics
    i_mo = args.i_mo
    ref_line_path = args.ref_line_path
    out_path = args.out_path
    

    '''
    Load configuration
    '''
    if seed == 0:
        CONFIG_FILE = "config_clean.yaml"
    else:
        CONFIG_FILE = "config_error.yaml"
    MASK_FILE = f"scripts/madx/mask_all_errors_{version}_{optics}.madx"

    with open(CONFIG_FILE) as f:
        config = yaml.safe_load(f)

    errors = config["errors"]
    sequence_name = "lhcb1" if beam == 1 else "lhcb2"

    print(f"Beam:     {beam} (mylhcbeam={'1' if beam == 1 else '4'})")
    print(f"Sequence: {sequence_name}")
    print(f"Seed:     {config['seed']}")
    print(f"Errors:")
    for err_type, enabled in errors.items():
        print(f"  {err_type:10s}: {'ON' if enabled else 'OFF'}")


    '''
    Generate error tables (only for enabled error types)
    '''
    arcs_errors    = errors["arcs"]
    D1_errors      = errors["D1"]
    D2_errors      = errors["D2"]
    MCBRD_errors   = errors.get("MCBRD", False)
    triplet_errors = errors["triplet"]

    if D1_errors:
        print("\n--- Generating D1 error tables ---")
        build_d1_tables(config, output_dir="tables/D1_tables")

    if D2_errors:
        print("\n--- Generating D2 error tables ---")
        build_d2_tables(config, output_dir="tables/D2_tables")

    if MCBRD_errors:
        print("\n--- Generating MCBRD error tables ---")
        build_mcbrd_tables(config, output_dir="tables/MCBRD_tables", seed=config["seed"])

    if triplet_errors:
        print("\n--- Generating triplet error tables ---")
        build_triplet_tables(config, output_dir="tables/IT_tables")

    print("\nTable generation complete.")


    '''
    Run MADX
    '''
    print("Running MAD-X mask...")
    madx = run_madx_mask(MASK_FILE, params={
        "mylhcbeam":         1 if beam == 1 else 4,
        "par_oct_current":   i_mo,
        "ver_hllhc_optics":  version,
        "par_myseed":        seed,
        "par_on_errors_LHC": int(arcs_errors),
        "on_D1":             int(D1_errors),
        "on_D2":             int(D2_errors),
        "on_MCBRD":          int(MCBRD_errors),
        "on_IT":             int(triplet_errors),
    }, verbose=True, acc_models_lhc=config["acc_models_lhc"])
    print("MAD-X execution complete.")


    '''
    Convert to xsuite line and save
    '''
    print(f"\nConverting sequence '{sequence_name}' to xtrack line...")
    line = madx_to_xtrack_line(madx, sequence_name)
    print(f"  {len(line.element_names)} elements")

    # Correct
    # if seed != 0:
    #     if ref_line_path == "None":
    #         raise ValueError("Reference line must be provided for spool piece " \
    #         "assignment and final correction!")
    #     line_ref = xt.load(ref_line_path)
    #     tw_ref = line_ref.twiss()
        
    #     assign_spool_pieces(line, beam)

    #     match_tune_chrom(line, tw_ref.qx, tw_ref.qy, tw_ref.dqx, tw_ref.dqy, tol=[1e-4, 5e-5, 1e-5, 5e-6])
    #     match_coupling(line, beam, c_minus=0.001, tol=1e-5)
    #     match_tune_chrom(line, tw_ref.qx, tw_ref.qy, tw_ref.dqx, tw_ref.dqy, tol=[1e-4, 5e-5, 1e-5, 5e-6])
    
    # Build output filename
    err_parts = [k for k, v in errors.items() if v]
    err_tag = "_".join(err_parts) if err_parts else "no_errors"
    OUTPUT_FILE = f"{out_path}/hllhc_b{beam}_v{version}_{optics}_imo{i_mo:.0f}_{err_tag}.json"

    line.to_json(OUTPUT_FILE)
    file_size = os.path.getsize(OUTPUT_FILE) / 1e6
    print(f"Saved to: {OUTPUT_FILE} ({file_size:.1f} MB)")


    '''
    Clean up
    '''
    print("Cleaning up temporary files...")
    cleanup()
    print("Done.")


if __name__ == "__main__":
    main()