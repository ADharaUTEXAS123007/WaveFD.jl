WaveFD.jl
=========
WaveFD.jl contains reference implementations of a collection of single time step wavefield propagation codes. The code is written in a mix of Julia and C++. For all propagators, we expose Julia and C++ APIs for propagating wavefields.  

# Propagator names
The propagator names are descriptive of physics and dimensionality, and described in the table below.

| Propagator | Dim  | Isotropy | Variable Density | Dissipation Only Attenuation | Time Order | Domain |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
Prop2DAcoIsoDenQ_DEO2_FDTD | 2D | Iso | Yes | Yes | 2 | Time |
Prop2DAcoVTIDenQ_DEO2_FDTD | 2D | VTI | Yes | Yes | 2 | Time |
Prop2DAcoTTIDenQ_DEO2_FDTD | 2D | TTI | Yes | Yes | 2 | Time |
Prop3DAcoIsoDenQ_DEO2_FDTD | 3D | Iso | Yes | Yes | 2 | Time |
Prop3DAcoVTIDenQ_DEO2_FDTD | 3D | VTI | Yes | Yes | 2 | Time |
Prop3DAcoTTIDenQ_DEO2_FDTD | 3D | TTI | Yes | Yes | 2 | Time |

# Micro Benchmark
You can run `benchmark/benchmarks.jl` to reproduce this micro benchmark on your hardware of choice. We use `PkgBenchmark.jl` to average over runs and generate repeatable timing information 

Please see `benchmark/README.md` for details. 

#### Architecture 
`Intel(R) Xeon(R) Platinum 8168 CPU @ 2.70GHz, 2 sockets, 44 cores`

Azure VM SKU: `Standard-HC44rs`

#### Model size
| Dimension | size X | size Y | size Z |
|:---------:|:------:|:------:|:------:|
| 2D        |   1000 |    --- |    500 |
| 3D        |   1000 |    300 |    500 |

#### Date of experiment 
October 6, 2020

#### Throughput (megacells/second)
|------| 1 threads | 2 threads | 4 threads | 8 threads | 16 threads | 32 threads | 44 threads|
|------| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 2DAcoIsoDenQ_DEO2_FDTD | 285.56 MC/s (3.35 %) | 547.28 MC/s (4.04 %) | 1037.57 MC/s (0.98 %) | 1847.97 MC/s (5.67 %) | 3020.84 MC/s (2.95 %) | 2593.60 MC/s (4.14 %) | 2006.65 MC/s (6.07 %)|
| 2DAcoVTIDenQ_DEO2_FDTD | 51.97 MC/s (1.36 %) | 102.39 MC/s (0.92 %) | 198.96 MC/s (1.91 %) | 376.68 MC/s (1.83 %) | 659.94 MC/s (1.21 %) | 709.06 MC/s (2.99 %) | 581.64 MC/s (5.33 %)|
| 2DAcoTTIDenQ_DEO2_FDTD | 35.46 MC/s (1.45 %) | 68.81 MC/s (1.56 %) | 133.71 MC/s (1.41 %) | 255.49 MC/s (1.84 %) | 472.84 MC/s (1.48 %) | 744.84 MC/s (2.61 %) | 691.24 MC/s (3.77 %)|
| 3DAcoIsoDenQ_DEO2_FDTD | 104.39 MC/s (0.59 %) | 202.05 MC/s (0.87 %) | 391.25 MC/s (0.57 %) | 710.36 MC/s (1.53 %) | 1000.61 MC/s (1.15 %) | 1625.00 MC/s (1.42 %) | 1961.90 MC/s (2.59 %)|
| 3DAcoVTIDenQ_DEO2_FDTD | 31.29 MC/s (0.66 %) | 61.29 MC/s (0.40 %) | 119.93 MC/s (0.34 %) | 226.95 MC/s (0.39 %) | 375.16 MC/s (1.14 %) | 692.11 MC/s (2.32 %) | 893.94 MC/s (3.40 %)|
| 3DAcoTTIDenQ_DEO2_FDTD | 26.95 MC/s (6.11 %) | 53.57 MC/s (0.26 %) | 105.18 MC/s (0.32 %) | 199.80 MC/s (1.78 %) | 339.85 MC/s (1.33 %) | 621.49 MC/s (1.75 %) | 800.58 MC/s (2.72 %)|

# See Also
* The `JetPackWaveFD.jl` package
* *Self-adjoint, energy-conserving second-order pseudoacoustic systems for VTI and TTI media for reverse time migration and full-waveform inversion* (2016)
https://library.seg.org/doi/10.1190/segam2016-13878451.1

