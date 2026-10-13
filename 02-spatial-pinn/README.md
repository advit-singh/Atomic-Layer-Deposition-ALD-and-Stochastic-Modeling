# Spatial network

`spatial_pinn_training.ipynb` trains one PhysicsNeMo `FullyConnected` network on thickness at 100 points along a trench period. Six layers, width 128. Inputs are temperature, pulse time, cycle count, and TMA mole fraction. The run uses 3,000 epochs, Adam at 2e-3 with weight decay 1e-5, and an 80/20 split.

The loss is data MSE plus 0.5 times three penalties: the derivative of top-fin thickness with respect to cycle count, matched to the Cantera growth per cycle; the second derivative of the mean thickness with respect to cycle count; and the thickness predicted when the cycle count is set to 0. The top-fin points are the outer 5% of the 100 indices at each end.

The full table, `ald_spatial_dataset_n2500.csv`, has 2,500 rows and is not in git. `dataset/sample_n2500_first20.csv` is the header and the first 20 rows. In `dataset/dataset_statistical_summary.csv`, `thickness_sidewall_nm` is 0.0303030303030311 nm for all 2,500 rows. The standard deviation on that column is numerical noise. The other thickness columns do vary. A level-set film that keeps one sidewall thickness for every process condition is what the September discussion was about. The note in `04-stochastic-plan` is the follow-up.

`dataset/spatial_pointwise_metrics.csv` is the validation split, one row per spatial point. The mean of those 100 MAE values is 4.2143034188484305e-05 nm. Per-point R² runs from 0.9997907876968384 to 0.9998862147331238.

Also in `dataset/`: the training history, the correlation matrix, and the top 5% error-tail rows. The figures are the heatmaps, the pointwise error profile, the error-tail histogram, and the four-panel diagnostic. Weight files are not in git.
