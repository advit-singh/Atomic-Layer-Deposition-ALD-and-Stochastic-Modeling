# Stochastic plan

There is no kinetic Monte Carlo script in this repository.

On 21 September I went through the level-set results with the group. The surrogate fits the thickness field. The field is an ideal film. In the 2,500-row summary, `thickness_sidewall_nm` does not change with temperature, pulse time, or cycle count. The level-set growth does not put roughness on the wall.

The model I want next has two Monte Carlo stages, and I have not written either one. Test-particle Monte Carlo would estimate the flux along the trench. Kinetic Monte Carlo would grow the film at those sites, then the pair would repeat for each cycle. ViennaPS was already used here in level-set mode. The Monte Carlo path is a different mode of that library. Cantera would still supply the surface chemistry.

Papers I was reading, by the titles on the files:

- A graph-theoretical kinetic Monte Carlo framework for on-lattice chemical kinetics
- A local superbasin kinetic Monte Carlo method
- Accurate acceleration of kinetic Monte Carlo simulations through the modification of rate constants
- Coarse-grained stochastic processes and kinetic Monte Carlo simulators for the diffusion of interacting particles
- First principles study of the atomic layer deposition of alumina by TMA-H2O process
- Kinetic Monte Carlo study of the atomic layer deposition of zinc oxide
- Sticking probabilities of H2O and Al(CH3)3 during atomic layer deposition of Al2O3 extracted from their impact on film conformality
- The graph transformation approach for first passage times of absorbing Markov chains
- ViennaPS, a flexible framework for semiconductor process simulation

The PDFs are not in git.
