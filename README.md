# Al2O3 ALD surrogates

This work is conducted as part of the Khan Lab at Georgia Tech. The process is thermal Al2O3 from trimethylaluminum and water, on a FinFET fin 42 nm tall. [Cantera](https://cantera.org/) solves the surface reactions. [ViennaPS](https://github.com/ViennaTools/ViennaPS) grows the film. A [PhysicsNeMo](https://github.com/NVIDIA/physicsnemo) fully connected network is trained on those runs.

The folders follow the order of the work.

`01-level-set-twin` builds five fins, runs the surface mechanism, grows a level-set film over 100 conditions, and trains a small surrogate.

`02-spatial-pinn` predicts thickness at 100 points along one trench period. The table behind it has 2,500 rows. Git has the first 20.

`03-factorial-ablation` is 16 networks: four loss formulations and four training sizes. The comparison is `results/ablation_matrix_results.csv`.

`04-stochastic-plan` has no solver. It is the Monte Carlo growth model I have not written, and the papers I was reading for it.

For model P3 trained on 8,000 samples, the ablation table gives a global spatial R² of 0.9999998807907104 and a mean absolute error of 4.144771082792431e-05 nm. `eval_duration_sec` on that row is 0.000914. The error is pointwise thickness against the level-set data. In the 2,500-row summary, the column `thickness_sidewall_nm` is 0.0303030303030311 nm on every row. Conformality R² at 20 and 100 samples is negative in the ablation table. P0 at 20 samples is -83076.09375. P3 at 8,000 samples is 0.9580895900726318.

The full spatial tables, the per-run meshes, and the weight files are not in git. Run a notebook from its own directory. `BASE_PATH` is that directory.

`slides/6-16-deposition-simulation-and-modeling.pptx` is the 16 June deck. It times four public codes. HackerFab is a valve script around hardware that wants an Arduino: 13,920 s on the hardware path, purge 3,000 ms, pulse counts 8, 12, and 4. DREAM-ALD on a 500 µm by 10 µm by 1,000 µm channel ran 500 cycles in 229.29 s, with a reported mean thickness of 0.162 nm. Machball on an aspect-ratio-100 via and trench finished in under 2 s. The PhysicsNeMo slide is a Fourier neural operator on a Tesla T4, 20 epochs and 100 samples, inference 1.52 to 2.22 ms. Those runs use their own geometries. The FinFET notebooks are the later folders.

The deck in `slides/first-and-second-simulation-presentation.pdf` is the spatial check and the 16-model matrix. Where the deck and a CSV disagree, the CSV is the one to use.

`slides/8-24-slides.pdf` is the 24 August deck. It names the 100 thickness indices: fin peaks at 0 to 5 and 94 to 99, sidewalls at 6 to 39 and 61 to 93, trench bottom at 40 to 60. It prints the P3 error as 0.000041 nm and R² 0.999988. Those are the rounded ablation-table figures. The deck also says a profile evaluation went from minutes to under 0.001 ms per sample. The minutes half of that sentence is only in the deck. `training_and_inference_timings.csv` records 0.00056 ms for one P3 sample trained on 8,000 rows.
