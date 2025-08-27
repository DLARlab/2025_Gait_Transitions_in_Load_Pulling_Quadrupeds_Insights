# Gait Transitions in Load-Pulling Quadrupeds: Insights from Sled Dogs and a Minimal SLIP Model

## Overview

This repository contains MATLAB code to simulate and analyze galloping gait transitions in sled dogs pulling loads, based on a minimal Spring-Loaded Inverted Pendulum (SLIP) model with hybrid dynamics. The code reproduces experimentally observed gait sequences, tugline forces, and stride-to-stride transitions by optimizing model parameters. It supports both single-gait replication and multi-stride transition analysis, enabling exploration of the biomechanical mechanisms behind locomotor multistability under load.


<p align="center">
  <img src="/Fig1_Overview.png" alt="Model" width="80%" height="80%">
</p>


Contents:
* **Section 1 — Gait Statics:** GUI to visualize gait timing statistics (footfall diagrams, duty factors, phase lags) from **animal experiments** and **SLIP simulations**.
* **Section 2 — Single-Stride Replication:** GUI to reproduce a **single periodic stride** with the hybrid SLIP–load model and compare to experimental averages.
* **Section 3 — Gait-Transition Replication:** GUI to replicate **stride-to-stride transitions** (e.g., TranR ↔ RotL) using limb-specific swing-stiffness modulation and to compare simulated tugline forces vs. measured data.
* **Stored_Functions:** Shared utilities for dynamics, hybrid events, optimization, plotting, and GUI components.

Features:
* **Hybrid SLIP-Load Model:** Planar quadrupedal SLIP model with load coupling via a unilateral spring (tugline).
* **Experimental Data Integration:** Matches footfall timings and tugline forces from high-speed sled dog experiments.
* **Transition Analysis:** Identifies limb-specific stiffness modulation strategies for stride-to-stride transitions.
* **Visualization:** Generates stance-phase diagrams, tugline force plots, and animations.

The source code is released under a [BSD 3-Clause license](LICENSE).

**Authors:** Jiayu Ding*, Benjamin Seleb*, Heather J. Huson, Saad Bhamla†, Zhenyu Gan†  
*These authors contributed equally  
†Corresponding authors  
**Affiliations:**
- College of Engineering and Computer Science, Syracuse University, USA
- Interdisciplinary Graduate Program in Quantitative Biosciences, Georgia Institute of Technology, USA
- Department of Animal Science, Cornell University, USA
- School of Chemical and Biomolecular Engineering, Georgia Institute of Technology, USA

<p align="center">
  <img src="/Fig2_Model.png" alt="Model" width="75%" height="75%">
</p>

## Publications

This work is available on **arXiv** (preprint).

If you use this work in an academic context, please cite the following publication:

> **J. Ding, B. Seleb, H. J. Huson, S. Bhamla, and Z. Gan**,  
> *"Gait Transitions in Load-Pulling Quadrupeds: Insights from Sled Dogs and a Minimal SLIP Model"*,  
> [arXiv:2507.14727](https://doi.org/10.48550/arXiv.2507.14727), 2025.

```bibtex
@misc{ding2025sleddoggaits,
      title={Gait Transitions in Load-Pulling Quadrupeds: Insights from Sled Dogs and a Minimal SLIP Model}, 
      author={Jiayu Ding and Benjamin Seleb and Heather J. Huson and Saad Bhamla and Zhenyu Gan},
      year={2025},
      eprint={2507.14727},
      archivePrefix={arXiv},
      primaryClass={cs.RO},
      doi={10.48550/arXiv.2507.14727}
}
```

The source code and datasets for the simulations and figures are available at:  
[https://github.com/DLARlab/2025_Gait_Transitions_in_Load_Pulling_Quadrupeds_Insights](https://github.com/DLARlab/2025_Gait_Transitions_in_Load_Pulling_Quadrupeds_Insights)

Details on the biologging devices used in the study, their design files, and supporting analysis scripts are available at:  
[https://github.com/bhamla-lab/DoggyLogger](https://github.com/bhamla-lab/DoggyLogger)

## Requirements

This code requires MATLAB R2019b or later.

## Usage

### Quick start (GUI-based)
All sections ship with a MATLAB GUI and include datasets from **animal experiments** and **model simulations**.

1. **Add repo to MATLAB path** (or `cd` into the repo root).
2. **Launch a section by function name** from the Command Window:
   - `Section1_Gait_Statics`
   - `Section2_Single_Stride_Replication`
   - `Section3_Gait_Transition_Replication`
3. In the GUI, click **Select Folder** and choose a dataset directory (experimental or simulation).  
4. Use the **data selection dropdown** to pick a specific dog/stride set or simulated run.  
5. Press the on-screen **Run/Update** controls to render plots (footfall bars, phase lags), tugline force traces, animations, and summary metrics.

### Notes
- Datasets are organized per section; both **experimental** and **simulation** examples are provided.
- Core utilities live in `Stored_Functions/` and are loaded automatically by the section GUIs.
- MATLAB R2019b+ is required.

