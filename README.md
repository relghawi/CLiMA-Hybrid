# CLiMA-Hybrid Repository

## Repository to run CLiMA-Hybrid model

### To run CLiMA independently:
1) ``` cd CLiMA-Hybrid/Land  ```
2) Initialize Julia REPL, and activate Julia env, and download needed packages.
3) ``` cd examples```
4) Edit relevant dirs in ``` example.jl ``` for input dataset  ```debug.nc``` and ```debug.jld2```, and specify output path.
5) Also edit paths to include files in ```CLiMA-Hybrid/Land/src/StomataModels/model/empirical.jl``` that refer to the Hybrid model files.
6) ``` include("example.jl") ``` for standard CLiMA run.

### To run hybrid model:
1) ``` cd CLiMA-Hybrid/src/Hybrid_model  ```
2) Activate Julia env.
3) Edit paths for inputs and outputs in ```hyb_model_ex.jl```.
4) ``` include("hyb_model_ex.jl") ``` for standard hybrid model run.

### To run hybrid model in CLiMA:
1) cd  ``` CLiMA-Hybrid/Land/examples ```
2) After making edits to input and output dirs, run  ``` include("example_hyb.jl") ```
