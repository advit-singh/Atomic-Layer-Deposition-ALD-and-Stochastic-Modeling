# Level-set twin

`finfet_level_set_twin.ipynb` is the first simulation. It builds the fin, writes a Cantera mechanism, grows a film with ViennaPS, and trains a small PhysicsNeMo network on the scalar outputs.

The geometry is five fins. Height 42 nm, base width 8 nm, top width 6 nm, pitch 42 nm, length along the fin 120 nm. The mesh in `geometry/post_photolith_finfet.vtp` is that surface. The grown-film meshes are not in git. The preview and the cross-section figures are.

The mechanism is TMA + OH(s) and then water, producing Al2O3 and methane. Site density is 1.66e-9 mol/cm². `chemistry/growth_per_cycle.csv` is the saturation sweep at a TMA mole fraction of 1e-3. Growth per cycle there is 0.2435379418717915 nm at 150 °C and 0.2704345077812945 nm at 300 °C.

The process sweep is a Latin hypercube of 100 conditions, seed 42. Temperature 150 to 300 °C, pulse time 0.05 to 2 s, cycles 4 to 20, TMA mole fraction 5e-4 to 5e-3. ViennaPS is called as a level-set process with sticking probabilities 0.003 for TMA and 0.01 for water. The per-run surfaces from that loop are not in git. `dataset/ald_high_fidelity_dataset.csv` and `dataset/design_of_experiments.csv` are.

The network is `FullyConnected` with 4 layers and width 64. The import that ran with this notebook reported PhysicsNeMo 2.1.1. `dataset/surrogate_vs_groundtruth.csv` has 12 rows. On total thickness, the median relative error is 0.646% and the maximum is 12.864%.

Figures:

- `figures/post_photolith_finfet_preview.png`
- `figures/ald_cycle_progression.png`
- `figures/before_after_ald.png`
- `figures/trench_conformality_cross_section.png`
- `figures/surrogate_training_curve.png`
- `figures/surrogate_benchmark.png`
