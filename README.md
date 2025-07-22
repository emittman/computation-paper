# Computational strategy for fitting Dirichlet Process Mixture Model for Gene Expression Experiments
- Second chapter of my dissertation [[https://dr.lib.iastate.edu/server/api/core/bitstreams/42b94ef8-a562-4776-80b9-f0f59e000452/content]]
- The study explores a flexible Bayesian nonparametric hierarchical model for gene expression data.
- High-throughput RNA sequencing enables simultaneous measurement of relative expression for tens of thousands of genes (Wang et al., 2009).
- Researchers aim to identify genes involved in various biological phenomena using gene expression data.
- Limited sample sizes in typical studies result in low statistical power for detecting differences between experimental groups.
- Hierarchical modeling aids in stabilizing estimation and reducing false detection rates through partial pooling of information without heavy reliance on tuning parameters.
- Some existing hierarchical methods have unrealistic assumptions that can reduce the efficiency of information pooling.
- The proposed approach avoids these assumptions and employs a Gibbs sampler on a GPU to make posterior inference computationally feasible.
