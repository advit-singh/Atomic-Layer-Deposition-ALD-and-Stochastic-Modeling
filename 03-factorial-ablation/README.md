# Factorial ablation

Two notebooks and one generator.

`generate_ald_dataset.py` builds the spatial tables. `N_SAMPLES` in the script is 10,000. It writes chemistry, topography, and dataset folders next to the script. The topography meshes and the full tables are not in git.

`dataset/sample_n10000_first20.csv` is the first 20 rows of a file with 10,000 rows. `dataset/sample_n100000_first20.csv` is the first 20 rows of the file named `ald_spatial_dataset_n100000.csv`. That file has 4,400 rows, not 100,000.

`exp1_holdout_baseline.ipynb` trains one network with data MSE only. The split is up to 8,000 training rows and 2,000 holdout rows, seed 42. Six layers, width 128, 2,500 epochs, patience 300. `results/exp1_corrected_train_vs_test_metrics.csv` is that holdout.

On the 2,000-row holdout, global spatial R² is 0.9999998211860657 and mean absolute error is 5.0165887e-05 nm. The angstrom column on that row is 0.0005016589. Maximum pointwise error is 0.00052411854 nm. Fin-peak, sidewall, trench-bottom, and corner MAEs are 5.5315733e-05, 5.0302184e-05, 4.711389e-05, and 5.3248557e-05 nm.

The deck in `slides/first-and-second-simulation-presentation.pdf` prints this holdout R² as 0.999988 and the MAE as 0.0005016589 nm. That MAE figure matches the angstrom column, not the nm column. The peak error and the regional MAEs on the deck match the nm columns.

`factorial_ablation.ipynb` trains 16 networks. Training sizes are 20, 100, 1,000, and 8,000. The four losses are:

- P0 is data MSE.
- P1 adds a penalty on thickness when the cycle count is set to 0.
- P2 adds the top-fin derivative matched to growth per cycle.
- P3 uses a weighted data MSE (sidewall indices 6 to 39 and 61 to 93 at 1.5, indices 40 to 60 at 3), the P1 and P2 terms, a second derivative of mean thickness with respect to cycles that ramps in over the first 10 cycles, and 0.1 times the square of the discrete second difference along the trench. The 24 August deck uses the same cuts and calls indices 0 to 5 and 94 to 99 the fin peaks.

The physics terms are multiplied by 0.5. The network is again six layers of width 128. Adam at 2e-3, weight decay 5e-5, gradient clip 1, up to 2,500 epochs, patience 300.

`results/ablation_matrix_results.csv` is the comparison. P3 at 8,000 samples: R² 0.9999998807907104, MAE 4.144771082792431e-05 nm, RMSE 5.6063665965114044e-05 nm, 95th percentile absolute error 0.0001075512554962188 nm, conformality R² 0.9580895900726318, `eval_duration_sec` 0.000914. P0 at the same size has MAE 4.725137114292011e-05 nm and conformality R² 0.9647344350814819. P0 at 20 samples has conformality R² -83076.09375. Global spatial R² is already 0.999915599822998 at that size, which is a poor description of conformality.

`results/training_and_inference_timings.csv` times the same 16 runs. P3 at 8,000 samples trained for 76.5402 s over 2,500 epochs on the machine in `results/system_hardware_environment_spec.json`: Python 3.12.13, PyTorch 2.11.0+cu128, one Tesla T4, 12.67 GB of RAM. Test-set inference time on that row is 0.00112 s. Single-sample latency in the same file is 0.00056 ms. Those are different columns from `eval_duration_sec`.

The master metric CSV and JSON repeat the 16-model table with more fields. `pointwise_spatial_profiles_all_16_models.csv` is the per-point profile. The figures are the four-panel ablation plot, the heatmaps, and the error-tail plots for the 2,000-row and 8,000-row checks. Weights are not in git.
