# Case study 3

## Getting started

### Installing Julia

See the official Julia documentation for instructions of how to install: https://docs.julialang.org/en/v1/manual/getting-started/#man-getting-started.

**(!!!) This case study requires you to have installed Julia 1.6.**

### Installing dependencies

First we clone the [`Epimap.jl`](https://github.com/epimap/Epimap.jl-public) repo using your favorite method. Using the shell:

``` shell
git clone https://github.com/epimap/Epimap.jl-public Epimap.jl
```

In particular we want the `#paper` branch, which for example can be achieve by executing the following command:

``` shell
git --git-dir=Epimap.jl/.git checkout paper
```

Now we can start the Julia REPL with the `Epimap.jl` project activated:

```shell
julia --project=Epimap.jl
```

Once in the repl is open, hit `]` to enter the REPL for Julia's package manager, and then type `instantiate` to install the dependencies.

``` julia
pkg> instantiate
```

## Getting the data

See [`epimap-data`](https://github.com/epimap/epimap-data) for how to obtain the data.

Once that has been done, set the environment variable `EPIMAP_DATA=/path/to/epimap-data/processed_data`. On Linux or MacOS, this can be done by

``` julia
export EPIMAP_DATA=/path/to/epimap-data/processed_data
```

before starting the repl using

``` julia
julia --project=Epimap.jl
```

## Inference

In this case study there are two models we'd run inference for:
1. `rmap`: observations are raw daily counts.
2. `rmap_debiased`: observations are debiased prevelance estimates obtained as described in the paper.

Note that both of these can take up to multiple days to run, depending on your hardware.
The experiments present in the paper took ~24hrs to run using the following setup:

``` julia
julia> versioninfo()
Julia Version 1.6.5
Commit 9058264a69 (2021-12-19 12:30 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Core(TM) i7-6850K CPU @ 3.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
```

The output of the inference scripts will be a folder of the format `Epimap.jl/intermediate/TIMESTAMP-BRANCH-COMMIT`, which contains everything related to the run.

*Note: that the inference results can vary between devices due to how the RNG is implemented.*

### `rmap`

```julia
julia> using DrWatson

julia> include(scriptsdir("run_rmap.jl"))
```

### `rmap_debiased`

```julia
julia> using DrWatson

julia> include(scriptsdir("run_debiased.jl"))
```

## Post-processing

Once inference has completed, we can compute several different statistics/"generated quantities" from the inference results. This is done using the `Epimap.jl/scripts/generate_outputs_*.jl*` scripts.

The resulting outputs can be found in a directory named `out` within the directory for the run, i.e. `Epimap.jl/intermediate/TIMESTAMP-BRANCH-COMMIT/out`. This is later used to produce the visualizations seen in the paper.

### `rmap`

``` julia
julia --project=Epimap.jl Epimap.jl/scripts/generate_outputs_rmap.jl Epimap.jl/intermediate/DIRECTORY-FOR-RMAP-RUN
```

### `rmap_debiased`

``` julia
julia --project=Epimap.jl Epimap.jl/scripts/generate_outputs_debiased.jl Epimap.jl/intermediate/DIRECTORY-FOR-RMAP_DEBIASED-RUN
```

## Figures

In the analysis, comparisons are made to an SIR model fitted to the debiased prevalence estimates from `case_study1`. This data can be found in the `data/` folder in this directory.

### Figure 8: geographic visualization

``` shell
julia --project=Epimap.jl scripts/mapviz.jl \
      --verbose \
      -drop-missing \
      --out=figures/figure-8.png \
      --date=2020-12-04 \
      --figure-size="(1500, 500)" \
      "Epimap=path/to/run-for-epimap-without-debiased/out/Rt.csv" \
      "Epimap debiased=path/to/run-for-epimap-with-debiased/out/Rt.csv" \
      "Debiased=path/to/output-for-case-study-1/IR_for_interop.csv"
```

![Figure 8](./figures/figure-8.png)

### Figure 9: Rt comparison between selected LTLAs

First we create a unified dataframe from the Rt quantities inferred for the different models:

``` shell
julia --project=Epimap.jl Epimap.jl/scripts/combine-with-case-study.jl --out=path/to/store/result path/to/debiased-Rt path/to/epimap-Rt path/to/epimap-debiased-Rt
```

Then we can run the script in `case_study1/06c_epimap_comparison_plots.R** on the resulting dataframe.

*Note: an example of such combined output can be found at `outputs/Rt-combined.csv`.*

### Figure 10: Internal vs. External "infection pressure"

``` shell
julia --project=Epimap.jl Epimap.jl/scripts/mapviz.jl \
      Epimap.jl/intermediate/DIRECTORY-FOR-RUN/out/Z_outside_portion.csv \
      --bounds="(0.0,1.0)" \
      --date=2020-12-04 \
      --column=Z_50 \
      --out=figures/figure-10.png \
      --drop-missing
```

![Figure 10](./figures/figure-10.png)
