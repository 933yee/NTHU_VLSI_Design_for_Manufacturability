# NTHU VLSI Design for Manufacturability

Coursework on algorithms that make advanced-node layouts easier and more reliable to manufacture. The assignments study lithography-aware optimization, redundant vias and fills, multiple-patterning decomposition, placement, and e-beam stencil planning; the final project implements a power-staple insertion optimizer in C++.

## Results

![HW1 grade](https://img.shields.io/badge/HW1-82-green)
![HW2 grade](https://img.shields.io/badge/HW2-98-green)
![HW3 grade](https://img.shields.io/badge/HW3-100-brightgreen)
![HW4 grade](https://img.shields.io/badge/HW4-79-green)
![Final grade](https://img.shields.io/badge/Final-88.14-green)

## Topics covered

| Work | Main topics |
| --- | --- |
| [HW1](HW1/) | Lithography limits and OPC, MANA wire-length matching, redundant-via conflict graphs, metal-fill optimization |
| [HW2](HW2/) | Minimum implant area repair, cell flipping/swapping graphs, ML hotspot features, double-patterning decomposition |
| [HW3](HW3/) | Triple-patterning decomposition, color-aware placement, stitch reduction, grid-based detailed routing |
| [HW4](HW4/) | Core-mask constraints in CNF/ILP, n-wise sampling, e-beam character projection, guiding-template optimization |
| [Final](Final/) | Triple-row dynamic programming for cell relocation and power-staple insertion |

Each homework directory contains the original assignment and the submitted solution report, including derivations, graph constructions, algorithm comparisons, and literature-based analysis.

## Final project — Power staple insertion

The final project implements the triple-row dynamic-programming method described by the course paper. Its objective is to maximize legal VDD/GND staple insertion while relocating standard cells within their displacement limits.

### Algorithm

1. Parse chip bounds, rows, sites, cell types/pins, initial placement, and displacement constraints.
2. Sort the cells in each row and build a site-to-cell map.
3. Process three rows at a time with a compact DP state containing the current site and placement progress for all three rows.
4. Enumerate deferral and single/multi-row placement transitions.
5. Reject transitions that cause overlap, excess displacement, pin blockage, anti-parallel conflict, staggering, alignment, or same-row violations.
6. Score each legal power-staple case and keep one state per compact key.
7. Backtrack from the highest-profit state, update cell positions, and continue with the next row group.
8. Balance VDD/GND counts and write the final cell/staple placement.

The implementation also includes an SVG debug renderer. For the largest cases, the submitted version reduces the effective displacement range to control DP state growth.

### Build and run

~~~bash
cd Final
g++ -O3 -std=c++17 main.cpp -o dfm_final
./dfm_final <input-file> <output-file>
~~~

The original benchmark inputs and checker are course material and are not included here. Input/output fields follow the format in [the project specification](Final/DFM%20project.pdf).

### Reported performance

The submitted report evaluates eight testcases, from 1,736 to 798,960 inserted staples, with measured runtimes of 16–520 seconds and peak memory below 1 GB in the reported environment. Full methodology and result figures are available in the [final report](Final/110062222_DFM_report.pdf).

## Repository structure

~~~text
HW1/ ... HW4/   Assignment specifications and solution reports
Final/main.cpp  C++ dynamic-programming implementation
Final/*.pdf     Project specification and implementation report
~~~

## Notes

- Only the final project contains executable source code; HW1–HW4 are theory/analysis assignments documented in PDF form.
- Reproducing final scores requires the course benchmark set, checker, and the same runtime limits.
