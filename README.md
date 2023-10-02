# CLiMA-Hybrid Repository

## Repository to run CLiMA-Hybrid model

### To run CLiMA independently:
1) ``` cd CLiMA-Hybrid/Land  ```
2) Initialize Julia REPL, and activate Julia env, and download needed packages.
3) ``` cd examples```
4)  ``` include("example.jl") ``` for standard CLiMA run.

### To run hybrid model:
1) ``` cd CLiMA-Hybrid/src/Hybrid_model  ```
2) Activate Julia env.
3) ``` include("hyb_model_ex.jl") ``` for standard hybrid model run.

### To run hybrid model in CLiMA:
1) cd  ``` CLiMA-Hybrid/Land/examples ```
2) Run  ``` include("example_hyb.jl") ```
