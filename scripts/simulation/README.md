# Simulation for Compound Het pipeline

For this section we use real genotypes and annotation with simulated phenotypes. We partition the heritability of the phenotype into three components, 1) a non-coding polygenic component, 2) a rare variant component and finally 3) a compound het component. We simulate under a combination of the infinitesimal model and the spike-slab model. We simulate combinations of phenotypes with the following option:

| Options | Description |
| --- | --- |
| h2 | We simulate the model under different heritabilities (0, 0.01, 0.05, 0.1) |
| pi | We simulate variying degree of polygenicity (0.01, 0.1, 1) |
| annotations | We simulate conditioned on underlying variant annotation (i.e. polygenicity for non-coding variants, spike-and-slab for coding variants and genes) |
| beta | We add fixed effects to gene effects (Probably something in the range of 0.001, 0.01, 0.01. |

