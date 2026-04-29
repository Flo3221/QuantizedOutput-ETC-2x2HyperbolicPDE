# QuantizedOutput-ETC-2x2HyperbolicPDE

## Description

QuantizedOutput-ETC-2x2HyperbolicPDE is a MATLAB/Simulink project developed to illustrate the concepts presented in the paper
**"Output Quantization Compensation in Event-Triggered Backstepping Control of 2×2 Hyperbolic Systems"** (submitted).

This project simulates a 2×2 first-order linear hyperbolic PDE system under **observer-based boundary control**, where the collocated output measurement is subject to **dynamic output quantization**. The control and observer gains are computed via the backstepping method. The control input is updated only at **event-triggered** (ETC) or **self-triggered** (STC) time instants, reducing actuation solicitation while guaranteeing semiglobal exponential stability in the sup-norm.

This project simulates:

- A 2×2 linear hyperbolic PDE system under **observer-based event-triggered boundary control** with dynamic output quantization (`q2x2sim.slx`).
- A 2×2 linear hyperbolic PDE system under **observer-based self-triggered boundary control** with dynamic output quantization (`selft_of_q2x2sim.slx`).

For detailed equations, backstepping kernel derivations, stability results, and parameter choices, please refer to the associated paper and ([mathematical details](https://tucgr-my.sharepoint.com/:b:/g/personal/fkoudohode_tuc_gr/IQAnmUwTMskQR7sLvugJRk1HAV8_9Rr44MHj5R2s6h0ikRI?e=99WzEn)).

## Requirements

To run this project, you will need:

- MATLAB R2023b or later (Symbolic Math Toolbox + Simulink required).

## Installation

Follow these steps to set up the project:

1. Download the project files from [QuantizedOutput-ETC-2x2HyperbolicPDE GitHub Repository](https://github.com/Flo3221/QuantizedOutput-ETC-2x2HyperbolicPDE).
2. Extract the contents to a directory of your choice.
3. Open MATLAB and navigate to the project directory using the `cd` command:

   ```
   cd /path/to/QuantizedOutput-ETC-2x2HyperbolicPDE
   ```

## Usage

### Order of execution

Run the scripts in the following order after downloading the project:

**Step 1 — Compute kernels and build the system matrices:**

```matlab
run('q2x2simu.m')
```

Defines system parameters, computes the backstepping control kernels `K` (`kernelsc.m`),
observer kernels `P` (`kernelso.m`), and their inverses `L` (`kernelscinv.m`),
`R` (`kernelsoinv.m`) via polynomial power series of order `N = 13`. Builds the
finite-difference matrices `A`, `Ao`, `B`, `Bo`, `C` and sets the initial conditions.

**Step 2 — Compute theoretical bounds:**

```matlab
run('bounds.m')
```

Evaluates all theoretical bounds (`a0`, `C1`, `Cth`, `Cp`, `Kq`, `M1`, `M2`,
`Cetmax`, `C0`), verifies feasibility conditions from the paper, and saves
`stc_constants.mat` for use inside the Simulink models.

**Step 3a — Run the ETC closed-loop simulation:**

```matlab
out = sim('q2x2sim.slx');
```

Runs the observer-based event-triggered closed-loop. Then generate figures:

```matlab
run('plots.m')
```

**Step 3b — Run the STC closed-loop simulation:**

```matlab
out = sim('selft_of_q2x2sim.slx');
```

Runs the observer-based self-triggered closed-loop. Then generate figures:

```matlab
run('STCplots.m')
```

**Step 4 — Inter-event time histogram (optional):**

```matlab
run('histogrametc.m')
```

Simulates 100 initial conditions over `[0, 15] s` and produces the density of
inter-event times `τ_{k+1} − τ_k` on a logarithmic scale, for two values of
`ν₀` and two values of `Ω`.

### Functions

QuantizedOutput-ETC-2x2HyperbolicPDE includes the following key functions.

#### Kernel computation

- `kernelsc.m`:
  Computes the backstepping **control kernels** K = (K^uu, K^uv, K^vu, K^vv) via a
  polynomial power series of order N, used to design the control gain vector **k**.
- `kernelso.m`:
  Computes the backstepping **observer kernels** P = (P^uu, P^uv, P^vu, P^vv) via a
  polynomial power series of order N, used to design the output injection gains
  **p₁**, **p₂**.
- `kernelscinv.m`:
  Computes the **inverse control kernels** L = (L^αα, L^αβ, L^βα, L^ββ), used to
  obtain the output-feedback control gain vector N_α, N_β.
- `kernelsoinv.m`:
  Computes the **inverse observer kernels** R, used in the theoretical bounds.

#### Setup and bounds

- `q2x2simu.m`:
  Defines system parameters (λ₁, λ₂, c₁, c₂, η, ρ), builds the finite-difference
  matrices A, Ao, B, Bo, C, computes all kernels (order N = 13), and sets initial
  conditions.
- `bounds.m`:
  Evaluates all theoretical bounds, verifies feasibility conditions from the paper,
  and saves `stc_constants.mat` for the Simulink models.

#### Simulation models

- `q2x2sim.slx`:
  Simulink model implementing the observer-based **event-triggered** closed-loop.
  Implements the plant ODE, observer ODE, dynamic quantizer on y(t) = u(t,1) with
  zoom ν(t) = ν₀ Ω^{⌊t/T⌋}, the ETC triggering mechanism, and the ZOH hold
  U_d(t) = U_nom(τ_k).
- `selft_of_q2x2sim.slx`:
  Same plant and observer, but the next trigger time is **pre-computed** at τ_k
  via the inter-execution time formula Γ(τ_k), so no continuous monitoring is needed.
- `q2x2simETCsup.slx`:
  Supplementary ETC Simulink model.

#### Plotting

- `plots.m`:
  Generates all ETC figures: state surface v(t,x), Lyapunov proxy Ψ(t), nominal
  vs. event-triggered control U_nom(t) and U_d(t), output y(t) vs. quantized output
  y_q(t) with zoom range ±Mν(t), and inter-execution time function r(τ_k).
- `STCplots.m`:
  Generates the same figures for the STC case, plus the inter-execution time r(τ_k)
  evaluated at each trigger instant τ_k as a stem plot.
- `histogrametc.m`:
  Produces the density of inter-event times τ_{k+1} − τ_k (log scale) over 100
  initial conditions, comparing two values of Ω and two values of ν₀.

## Examples

Refer to the following files for examples of how to use QuantizedOutput-ETC-2x2HyperbolicPDE:

- `q2x2sim.slx`
  (Observer-based event-triggered control with dynamic output quantization.)
- `selft_of_q2x2sim.slx`
  (Observer-based self-triggered control with dynamic output quantization.)

## Contributing

To contribute to QuantizedOutput-ETC-2x2HyperbolicPDE, please follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch for your feature or fix.
3. Make your changes and commit them.
4. Submit a pull request with a detailed description of your changes.

## License

This project is licensed under the CC BY-NC-ND license
([`LICENSE`](https://creativecommons.org/licenses/by-nc-nd/4.0/)).

## Contact

For questions or feedback, please contact [fkoudohode@tuc.gr](mailto:fkoudohode@tuc.gr).

# Acknowledgements

Funded by the European Union (ERC, C-NORA, 101088147). Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

## Cite this work

If you use this code in your research, please cite:

```
@article{koudohode2025quant2x2,
  title   = {Output Quantization Compensation in Event-Triggered Backstepping
             Control of {2$\times$2} Hyperbolic Systems},
  author  = {Koudohode, F. and Espitia, N. and Humaloja, J.-P. and Bekiaris-Liberis, N.},
  journal = {under review},
  year    = {2025}
}
```
