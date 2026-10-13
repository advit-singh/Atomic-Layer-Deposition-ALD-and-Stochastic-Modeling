import os
import gc
import time
import json
import traceback
import numpy as np
import pandas as pd
import pyvista as pv
import cantera as ct
import viennaps as vps
from scipy.stats import qmc
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from concurrent.futures import ProcessPoolExecutor, as_completed

N_SAMPLES = 10000
BASE_SAVE_PATH = os.path.dirname(os.path.abspath(__file__))
MECHANISM_YAML = os.path.join(BASE_SAVE_PATH, "chemistry", "ald_al2o3_mechanism.yaml")
DATASET_CSV = os.path.join(BASE_SAVE_PATH, "dataset", f"ald_spatial_dataset_n{N_SAMPLES}.csv")

BATCH_FLUSH_SIZE = 50
WORKER_MAX_TASKS = 100

FIN_HEIGHT_NM = 42.0
FIN_WIDTH_NM = 8.0
FIN_PITCH_NM = 42.0
SITE_DENSITY_MOL_CM2 = 1.66e-9
M_AL2O3_G_MOL = 101.96
RHO_AL2O3_G_CM3 = 3.1

def langmuir(t, theta_m, tau):
    return theta_m * (1.0 - np.exp(-t / tau))

def setup_environment():

    for folder in ["chemistry", "topography", "dataset"]:
        os.makedirs(os.path.join(BASE_SAVE_PATH, folder), exist_ok=True)

    mechanism_text = """
generator: ald-al2o3-mechanism-builder
cantera-version: 3.2.0
units: {length: cm, time: s, quantity: mol, activation-energy: kJ/mol}
phases:
  - name: gas
    thermo: ideal-gas
    elements: [Al, C, H, O, N]
    species: [TMA, H2O, CH4, N2]
    state: {T: 300.0, P: 1.0 atm}
  - name: Al2O3_bulk
    thermo: fixed-stoichiometry
    elements: [Al, O]
    species: [Al2O3s]
  - name: Al2O3_surface
    thermo: ideal-surface
    elements: [Al, C, H, O, N]
    species: [OH_s, O_Al_CH3_2_s]
    kinetics: surface
    reactions: all
    site-density: 1.66e-9
    adjacent-phases: [gas, Al2O3_bulk]
    state: {T: 300.0, coverages: {OH_s: 1.0, O_Al_CH3_2_s: 0.0}}
species:
  - name: TMA
    composition: {Al: 1, C: 3, H: 9}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [4.0, 0.02, -1.0e-5, 2.0e-9, -1.5e-13, -8000.0, 10.0]
        - [4.0, 0.02, -1.0e-5, 2.0e-9, -1.5e-13, -8000.0, 10.0]
  - name: H2O
    composition: {H: 2, O: 1}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [4.198, -2.036e-3, 6.520e-6, -5.487e-9, 1.771e-12, -30293.7, -0.849]
        - [3.033, 2.176e-3, -1.640e-7, -9.704e-11, 1.682e-14, -30004.3, 4.966]
  - name: CH4
    composition: {C: 1, H: 4}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [5.149, -0.01367, 4.918e-5, -4.847e-8, 1.666e-11, -10246.6, -4.643]
        - [1.911, 9.603e-3, -3.567e-6, 6.075e-10, -3.923e-14, -10099.2, 8.483]
  - name: N2
    composition: {N: 2}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [3.531, -1.237e-4, -5.030e-7, 2.435e-9, -1.409e-12, -1046.98, 2.968]
        - [2.953, 1.397e-3, -4.927e-7, 7.870e-11, -4.608e-15, -923.949, 5.872]
  - name: OH_s
    composition: {O: 1, H: 1}
    thermo:
      model: constant-cp
      T0: 300.0
      h0: -25.0 kJ/mol
      s0: 0.0 J/mol/K
      cp0: 0.0 J/mol/K
  - name: O_Al_CH3_2_s
    composition: {Al: 1, O: 1, C: 2, H: 6}
    thermo:
      model: constant-cp
      T0: 300.0
      h0: -120.0 kJ/mol
      s0: 0.0 J/mol/K
      cp0: 0.0 J/mol/K
  - name: Al2O3s
    composition: {Al: 2, O: 3}
    thermo:
      model: constant-cp
      T0: 300.0
      h0: -1675.7 kJ/mol
      s0: 50.9 J/mol/K
      cp0: 79.0 J/mol/K
    equation-of-state:
      model: constant-volume
      density: 3.1 g/cm^3
reactions:
  - equation: TMA + OH_s => O_Al_CH3_2_s + CH4
    rate-constant: {A: 3.6e+7, b: 0.5, Ea: 12.0}
  - equation: 2 O_Al_CH3_2_s + 3 H2O => 2 OH_s + Al2O3s + 4 CH4
    rate-constant: {A: 1.1e+8, b: 0.5, Ea: 18.0}
"""
    with open(MECHANISM_YAML, "w") as f:
        f.write(mechanism_text)

def compute_initial_z(x_coords):

    z_init = np.zeros_like(x_coords)
    hw = FIN_WIDTH_NM / 2.0
    for i, x in enumerate(x_coords):
        if abs(x) <= hw:
            z_init[i] = FIN_HEIGHT_NM
        else:
            z_init[i] = 0.0
    return z_init

def run_single_simulation(row_dict):

    pv.OFF_SCREEN = True
    if hasattr(ct, 'suppress_thermo_warnings'):
        ct.suppress_thermo_warnings()

    run_idx = int(row_dict['run_id'])

    try:

        surf_local = ct.Interface(MECHANISM_YAML, "Al2O3_surface")
        gas_local = surf_local.adjacent["gas"]
        T_K = row_dict['temperature_C'] + 273.15

        times = np.logspace(-2, 0.5, 10)
        coverages = []
        for t in times:
            surf_local.TP = T_K, ct.one_atm
            surf_local.coverages = {"OH_s": 1.0, "O_Al_CH3_2_s": 0.0}
            gas_local.TPX = T_K, ct.one_atm, {"TMA": row_dict['tma_mole_fraction'], "N2": 1.0 - row_dict['tma_mole_fraction']}
            surf_local.advance_coverages(t)
            coverages.append(surf_local["O_Al_CH3_2_s"].coverages[0])

        try:
            popt, _ = curve_fit(langmuir, times, coverages, p0=[1.0, 0.1], maxfev=2000)
            _, tau_tma_s = popt
        except Exception:
            tau_tma_s = 0.1

        surf_local.TP = T_K, ct.one_atm
        surf_local.coverages = {"OH_s": 1.0, "O_Al_CH3_2_s": 0.0}
        gas_local.TPX = T_K, ct.one_atm, {"TMA": row_dict['tma_mole_fraction'], "N2": 1.0 - row_dict['tma_mole_fraction']}
        surf_local.advance_coverages(row_dict['pulse_time_s'])

        theta_sat = surf_local["O_Al_CH3_2_s"].coverages[0]
        areal_mol_al2o3_cm2 = 0.5 * theta_sat * SITE_DENSITY_MOL_CM2
        gpc_nm = (areal_mol_al2o3_cm2 * M_AL2O3_G_MOL / RHO_AL2O3_G_CM3) * 1.0e7

        vps.setDimension(2)
        vps.setNumThreads(1)

        local_domain = vps.Domain()
        vps.MakeFin(
            domain=local_domain,
            gridDelta=1.0,
            xExtent=FIN_PITCH_NM,
            yExtent=FIN_HEIGHT_NM + 20.0,
            finWidth=FIN_WIDTH_NM,
            finHeight=FIN_HEIGHT_NM,
            baseHeight=0.0,
            periodicBoundary=True,
            makeMask=False,
            material=vps.Material.Si,
        ).apply()

        depo_model = vps.MultiParticleProcess()
        depo_model.addNeutralParticle(stickingProbability=0.003, label="TMA_flux")
        depo_model.addNeutralParticle(stickingProbability=0.01, label="H2O_flux")

        state = {"coverage": 0.0}
        def rate_fn(fluxes, material):
            tma_flux = fluxes[0]
            theta = state["coverage"]
            d_theta = (tma_flux * theta_sat - theta) / max(tau_tma_s, 1e-6) * 1e-3
            theta = min(max(theta + d_theta, 0.0), theta_sat)
            state["coverage"] = theta
            return gpc_nm * (d_theta / max(theta_sat, 1e-6))

        depo_model.setRateFunction(rate_fn)

        process = vps.Process(local_domain, depo_model)
        pulse_duration = tau_tma_s * 4.0
        process.setProcessDuration(pulse_duration)

        for _ in range(int(row_dict['num_cycles'])):
            state["coverage"] = 0.0
            process.apply()

        out_vtp = os.path.join(BASE_SAVE_PATH, "topography", f"run_{run_idx:04d}.vtp")
        local_domain.saveSurfaceMesh(out_vtp, True)

        final_mesh = pv.read(out_vtp)
        pts = final_mesh.points

        mask = (pts[:, 0] >= -FIN_PITCH_NM / 2.0) & (pts[:, 0] <= FIN_PITCH_NM / 2.0)
        x_raw = pts[mask, 0]
        z_raw = pts[mask, 1]

        sort_idx = np.argsort(x_raw)
        x_sorted = x_raw[sort_idx]
        z_sorted = z_raw[sort_idx]

        x_unique, unique_indices = np.unique(x_sorted, return_index=True)
        z_unique = z_sorted[unique_indices]

        f_interp = interp1d(x_unique, z_unique, kind='linear', bounds_error=False, fill_value="extrapolate")
        x_eval = np.linspace(-FIN_PITCH_NM / 2.0, FIN_PITCH_NM / 2.0, 100)
        z_final = f_interp(x_eval)
        z_initial = compute_initial_z(x_eval)

        vertical_thickness = z_final - z_initial

        trench_mask = (x_eval > 10.0) | (x_eval < -10.0)
        top_mask = (abs(x_eval) < 2.0)

        z_bottom = np.min(z_final[trench_mask]) if any(trench_mask) else 0.0
        z_top = np.max(z_final[top_mask]) if any(top_mask) else FIN_HEIGHT_NM

        thickness_bottom_nm = z_bottom - 0.0
        thickness_top_nm = z_top - FIN_HEIGHT_NM

        hw_mid = FIN_WIDTH_NM / 2.0
        mid_height_idx = np.argmin(np.abs(z_final - (FIN_HEIGHT_NM / 2.0)))
        x_sidewall = abs(x_eval[mid_height_idx])
        thickness_sidewall_nm = x_sidewall - hw_mid

        result = {
            "run_id": run_idx,
            "temperature_C": row_dict['temperature_C'],
            "pulse_time_s": row_dict['pulse_time_s'],
            "num_cycles": row_dict['num_cycles'],
            "tma_mole_fraction": row_dict['tma_mole_fraction'],
            "gpc_nm_per_cycle": gpc_nm,
            "theta_saturated": theta_sat,
            "thickness_top_nm": float(thickness_top_nm),
            "thickness_sidewall_nm": float(thickness_sidewall_nm),
            "thickness_bottom_nm": float(thickness_bottom_nm),
            "conformality_bottom_to_top": float(thickness_bottom_nm / max(thickness_top_nm, 1e-9)),
            "raw_x_coordinates": json.dumps(np.round(x_eval, 5).tolist()),
            "raw_z_coordinates": json.dumps(np.round(z_final, 5).tolist()),
        }

        for i in range(100):
            result[f"thickness_pt_{i:03d}"] = float(vertical_thickness[i])

        return result

    except Exception as e:
        print(f"[FAIL] Run {run_idx:04d} failed: {e}")
        traceback.print_exc()
        return None
    finally:
        pv.close_all()
        gc.collect()

def main():
    print(f"Initializing Parallel ALD Production Pipeline ({N_SAMPLES} Samples)")
    setup_environment()

    PARAM_BOUNDS = {
        "temperature_C": (150.0, 300.0),
        "pulse_time_s": (0.05, 2.0),
        "num_cycles": (4, 20),
        "tma_mole_fraction": (5e-4, 5e-3),
    }

    sampler = qmc.LatinHypercube(d=4, seed=42)
    lhs_unit = sampler.random(n=N_SAMPLES)
    bounds_lo = np.array([b[0] for b in PARAM_BOUNDS.values()])
    bounds_hi = np.array([b[1] for b in PARAM_BOUNDS.values()])
    lhs_scaled = qmc.scale(lhs_unit, bounds_lo, bounds_hi)

    design_df = pd.DataFrame(lhs_scaled, columns=list(PARAM_BOUNDS.keys()))
    design_df["num_cycles"] = design_df["num_cycles"].round().astype(int)
    design_df.insert(0, 'run_id', range(N_SAMPLES))

    completed_runs = set()
    write_header = not os.path.exists(DATASET_CSV)

    if os.path.exists(DATASET_CSV):
        try:
            existing_df = pd.read_csv(DATASET_CSV, usecols=['run_id'])
            completed_runs = set(existing_df['run_id'].values)
            print(f"[CHECKPOINT] Found {len(completed_runs)} completed runs. Resuming...")
        except Exception as e:
            print(f"[WARNING] Could not read checkpoint CSV: {e}. Resuming from empty set.")

    pending_rows = [row.to_dict() for _, row in design_df.iterrows() if int(row['run_id']) not in completed_runs]
    print(f"Remaining runs to execute: {len(pending_rows)} / {N_SAMPLES}")

    if not pending_rows:
        print("All runs are already completed.")
        return

    num_workers = max(1, os.cpu_count() - 1)
    print(f"Executing with {num_workers} parallel CPU workers...")

    t_start = time.time()
    results_buffer = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(run_single_simulation, row): row['run_id'] for row in pending_rows}

        for idx, future in enumerate(as_completed(futures), 1):
            res = future.result()
            if res is not None:
                results_buffer.append(res)

            if len(results_buffer) >= BATCH_FLUSH_SIZE or idx == len(pending_rows):
                df_batch = pd.DataFrame(results_buffer)

                written = False
                for attempt in range(5):
                    try:
                        df_batch.to_csv(DATASET_CSV, mode='a', header=write_header, index=False)
                        write_header = False
                        written = True
                        break
                    except PermissionError:
                        time.sleep(1.0)

                if not written:
                    print(f"[ERROR] Failed to write batch starting at run {results_buffer[0]['run_id']} due to Google Drive lock.")

                results_buffer.clear()
                print(f"[PROGRESS] Completed {idx}/{len(pending_rows)} pending runs ({(idx/len(pending_rows))*100:.1f}%).")

    t_end = time.time()
    print("\n" + "="*70)
    print(f"EXECUTION COMPLETED in {(t_end - t_start) / 3600:.2f} hours.")
    print(f"Data saved to {DATASET_CSV}")
    print("="*70)

if __name__ == '__main__':
    main()
