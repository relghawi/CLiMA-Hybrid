# CLiMA-Hybrid Repository

## Repository to run CLiMA-Hybrid model

Just activate the environment in the root folder [CLiMA-Hybrid] and do:

```julia
] add https://github.com/lazarusA/Land.git#single_site
```

Download these files

```
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/debug.jld2
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/debug.nc
```

into `ClimaLand_examples`.

### To run CLiMA example.jl:
0) Simply go to ClimaLand_examples/example.jl and start running line by line (or all the code) in vscode.
1) Or, initialize the Julia REPL, then activate and instantiate the Julia env and then do:
2) ``` cd ClimaLand_examples```
3)  ``` include("example.jl") ``` for standard REPL CLiMA run.

### To run hybrid model:
1) ``` cd Experiments  ```
2) ``` include("hyb_model_ex1.jl") ``` for standard hybrid model run.

### To run hybrid model in CLiMA:
1) cd  ``` CLiMA-Hybrid/Land/examples ```
2) Run  ``` include("example_hyb.jl") ```
