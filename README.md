# Process Engineering Toolkit (PET)

## Project Description

A lightweight modular toolkit for first-principles process design, equipment sizing, and operating window exploration of chemical unit operations.

PET aims to provide simple, physics-based modeling tools for rapid conceptual process development.

## Table of Contents

- [Motivation](#Motivation)
- [Project Status](#project-status)
    - [Planned Progression](#planned-progression)
- [Repository Structure](#repository-structure)

## Motivation

Modern process engineering tools are often:

* Heavy
* Black-box
* Difficult to extend

PET is designed as:

* Modular
* Explicit
* First-principles focused
* Easy to extend with custom unit operations

## Project Status

A general framework / system skeleton is now complete (up for revision at a later stage) and the first reactor has been added (Plug Flow Reactor - PFR). PFR has completed the test case of solving equilibrium conversion and new functionalities will added according to the Planned Progression below.

### Planned Progression

- [ ] Add new functionalities to PFR
    - [X] Determine required volume given target conversion
    - [X] Determine conversion profile
    - [X] Determine conversion for specific volume
- [ ] Add new general reactor functionalities
    - [ ] Automatic unit handling/conversion using pint
    - [ ] Additional rate laws
        - [ ] Custom
        - [ ] Michaelis-Menten
    - [ ] Plotting (including target folder for saving plots)
- [ ] Add additional units
    - [ ] CSTR
    - [ ] PBR
    - [ ] Batch
- [ ] Add recirculation
- [ ] Add serial and parallel reactor unit processing
- [ ] Add heating / cooling capabilities
- [ ] Add separation/filtration processes
- [ ] Add high-level infrastructure to manage parameter exploration
    - [ ] Grid search
    - [ ] Design of Experiments
    - [ ] Bayesian Optimization
- [ ] Add piping and pumps
- [ ] Add catalyst models


## Repository structure

```process_engineering_toolkit/
│
├── docs/                    # Documentation
├── notebooks/               # Example explorations and tutorials
├── plots/                   # Saved plots (auto-created by toolkit)
├── src/
│   └── process_engineering_toolkit/
│       ├── __init__.py
│       ├── stream.py        # Stream class definition
│       ├── reactions.py     # Reaction and Reactions classes
│       ├── rate_laws.py     # Rate laws: Only MassAction currently
│       └── reactors/
│           ├── __init__.py
│           ├── base.py      # Reactor base class
│           ├── pfr.py       # Plug Flow Reactor
│           ├── cstr.py      # CSTR (planned)
│           ├── batch.py     # Batch Reactor (planned)
│           └── pbr.py       # Packed Bed Reactor (planned)
│
├── README.md                # This file
└── pyproject.toml           # Project configuration and dependencies
```

<!-- ## Illustrations

Here you can include visuals like:

Example PFR concentration vs. conversion plots

Reactor network diagrams

High-level schematic of the PET workflow

Tip: Include images in docs/ or plots/ and reference them with Markdown:

![PFR Example](plots/pfr_example.png)

## Getting Started

# Clone the repo
git clone <repo-url>
cd process_engineering_toolkit

# Install in editable mode
pip install -e src/

# Run example notebooks or scripts
jupyter notebook notebooks/exploration.ipynb

## Dependencies

Python ≥ 3.11

NumPy

SciPy

Pint

Matplotlib (for plotting)

(Add any additional dependencies as they are added.)

## License -->