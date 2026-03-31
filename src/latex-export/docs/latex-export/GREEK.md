# Greek Alphabet — LaTeX Mapping Reference

This document defines the complete mapping from identifier names to LaTeX Greek letter commands.

## Matching rules

1. **Greedy matching by length:** always try longer names first to avoid partial matches (e.g., `varepsilon` must not match as `var` + `epsilon`).
2. **Exact match:** identifier is exactly a Greek letter name → replace entirely.
3. **Prefix + digit suffix:** identifier starts with a Greek letter name and the remainder is all digits → Greek letter with subscript. Example: `alpha1` → `\alpha_{1}`.
4. **No letter suffix:** if the remainder after the Greek prefix starts with a letter, do NOT treat as Greek. Example: `muM` → `\mathrm{muM}` (not `\mu M`). This prevents false positives like `delta` matching inside `deltaTime`.
5. **Case-sensitive:** `alpha` and `Alpha` are different entries.

## Lowercase Greek letters

Sorted by key length descending (for greedy matching):

| Key (length) | LaTeX | Unicode |
|---|---|---|
| `varepsilon` (10) | `\varepsilon` | ε |
| `vartheta` (8) | `\vartheta` | ϑ |
| `varsigma` (8) | `\varsigma` | ς |
| `upsilon` (7) | `\upsilon` | υ |
| `epsilon` (7) | `\epsilon` | ε |
| `omicron` (7) | `o` | ο |
| `varphi` (6) | `\varphi` | φ |
| `varrho` (6) | `\varrho` | ρ |
| `lambda` (6) | `\lambda` | λ |
| `sigmа` (5) | `\sigma` | σ |
| `omega` (5) | `\omega` | ω |
| `theta` (5) | `\theta` | θ |
| `gamma` (5) | `\gamma` | γ |
| `delta` (5) | `\delta` | δ |
| `kappa` (5) | `\kappa` | κ |
| `alpha` (5) | `\alpha` | α |
| `beta` (4) | `\beta` | β |
| `zeta` (4) | `\zeta` | ζ |
| `iota` (4) | `\iota` | ι |
| `eta` (3) | `\eta` | η |
| `tau` (3) | `\tau` | τ |
| `phi` (3) | `\phi` | φ |
| `chi` (3) | `\chi` | χ |
| `psi` (3) | `\psi` | ψ |
| `mu` (2) | `\mu` | μ |
| `nu` (2) | `\nu` | ν |
| `xi` (2) | `\xi` | ξ |
| `pi` (2) | `\pi` | π |

## Uppercase Greek letters

Only letters that differ visually from Latin have dedicated LaTeX commands. Others use `\mathrm{}`.

| Key | LaTeX | Note |
|---|---|---|
| `Gamma` | `\Gamma` | Γ |
| `Delta` | `\Delta` | Δ |
| `Theta` | `\Theta` | Θ |
| `Lambda` | `\Lambda` | Λ |
| `Xi` | `\Xi` | Ξ |
| `Pi` | `\Pi` | Π |
| `Sigma` | `\Sigma` | Σ |
| `Upsilon` | `\Upsilon` | Υ |
| `Phi` | `\Phi` | Φ |
| `Psi` | `\Psi` | Ψ |
| `Omega` | `\Omega` | Ω |
| `Alpha` | `\mathrm{A}` | Same as Latin A |
| `Beta` | `\mathrm{B}` | Same as Latin B |
| `Epsilon` | `\mathrm{E}` | Same as Latin E |
| `Zeta` | `\mathrm{Z}` | Same as Latin Z |
| `Eta` | `\mathrm{H}` | Same as Latin H |
| `Iota` | `\mathrm{I}` | Same as Latin I |
| `Kappa` | `\mathrm{K}` | Same as Latin K |
| `Mu` | `\mathrm{M}` | Same as Latin M |
| `Nu` | `\mathrm{N}` | Same as Latin N |
| `Omicron` | `\mathrm{O}` | Same as Latin O |
| `Rho` | `\mathrm{P}` | Same as Latin P (!) |
| `Tau` | `\mathrm{T}` | Same as Latin T |
| `Chi` | `\mathrm{X}` | Same as Latin X |

## Special constants

| Key | LaTeX |
|---|---|
| `PI` | `\pi` |
| `Inf` | `\infty` |
| `inf` | `\infty` |

## Variant forms

Some letters have variant forms commonly used in specific mathematical contexts:

| Standard | Variant | Common usage |
|---|---|---|
| `\epsilon` | `\varepsilon` | More common in analysis |
| `\theta` | `\vartheta` | Less common |
| `\rho` | `\varrho` | Less common |
| `\sigma` | `\varsigma` | Final sigma |
| `\phi` | `\varphi` | More common in physics |
