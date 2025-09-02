# Elevated plasma Tau-PT217 linked to decreased hippocampal functional connectivity in patients with knee osteoarthritis

🧠 Osteoarthritis, Plasma Tau, and Brain Connectivity: fMRI & Biomarker Analysis
📌 Background
Osteoarthritis (OA) has been increasingly associated with a higher risk of dementia, though the underlying biological mechanisms remain poorly understood. Recent work suggests that blood phosphorylated tau proteins (p-Tau), particularly Tau-PT217, may serve as early biomarkers for cognitive decline and Alzheimer’s disease.
In this project, we investigate:
- Plasma phosphorylated tau protein levels (Tau-PT217 and Tau-PT181)
- Resting-state fMRI functional connectivity, particularly hippocampal networks
- Cognitive performance in people with knee OA vs. age- and gender-matched pain-free controls

🧪 Data & Preprocessing
 - Plasma Biomarker Analysis
- Compared biomarker distributions between OA patients and controls using independent samples t-tests and ANCOVA (adjusting for age, sex, and education).
- fMRI Data Processing
- Structural MRI preprocessing: T1-weighted images processed with FreeSurfer for cortical reconstruction and volumetric segmentation.
- Functional MRI preprocessing: Motion correction, slice-timing correction, and spatial normalization (MNI space), Noise removal using CompCor and head motion regressors, Band-pass filtering to retain 0.01–0.1 Hz fluctuations, Spatial smoothing with a 6 mm FWHM Gaussian kernel.
- Seed-based functional connectivity: bilateral hippocampus seeds extracted.
- Statistical inference: Two-sample t-tests (OA vs. control) on seed-to-voxel maps, Cluster-level correction with FDR (q < 0.05).
- Cognitive Data: Neuropsychological tests included measures of memory, attention, and executive function,
- Correlations assessed between Tau-PT217 levels, hippocampal connectivity strength, and cognitive performance.

📊 Results
Blood biomarkers:
- OA patients had significantly higher plasma Tau-PT217 (p < 0.05). No group difference for Tau-PT181.
Functional Connectivity:
- OA patients showed reduced hippocampus–MCC connectivity (FDR-corrected p < 0.05).
- Connectivity strength negatively correlated with plasma Tau-PT217 levels.
Cognition:
- Elevated Tau-PT217 levels correlated with lower memory and executive function scores.

🚀 Key Takeaways
Plasma Tau-PT217 is a sensitive marker linked to both OA and brain connectivity alterations.
Reduced hippocampal connectivity with the MCC may represent an early mechanism of dementia risk in OA.
Findings support a biological pathway between knee osteoarthritis and neurodegeneration, opening new directions for early intervention.

Published paper can be found at https://doi.org/10.1016/j.brainres.2025.149478
