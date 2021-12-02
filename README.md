# SensitivityCompositeBinary

The file compositeBinary.R implements the methodology in Fogarty et. al (2017), entitled ``Randomization Inference and Sensitivity Analysis for Composite Null Hypotheses With Binary Outcomes in Matched Observational Studies,'' precisely as originally described. 

The file compositeBinaryCenter.R uses the same ideas as those presented in Fogarty et. al (2017) for conducting a sensitivity analysis with binary outcomes after matching, but uses a modified form for the test statistic inspired by results in Fogarty (2021+), ``Testing weak nulls in matched observational studies.'' Instead of looking at weighted differences in means, the code instead looks at the weighted difference in means minus a centering term. The resulting test statistic has worst-case expectation equal to zero under Fisher's sharp null, but markedly improves performance in sensitivity analysis even with heterogeneous effects relative to the original method implemented in compositeBinary.R. We recommend proceeding with this modification. 
