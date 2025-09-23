# ADMolarAvgVelocityMaterial

!syntax description /Materials/ADMolarAvgVelocityMaterial

## Overview

`ADMolarAvgVelocityMaterial` computes the molar-average velocity of a multicomponent mixture by
combining a provided bulk velocity with diffusive fluxes expressed in an \((N-1)\)-species basis.
The object assumes the species fluxes satisfy the closure constraint
\(\sum_i \mathbf{J}_i = 0\) and eliminates the final species flux when forming the molar-average
velocity
\(\mathbf{u}_\mathrm{m} = \mathbf{v} + \rho^{-1} \sum_{i=1}^{N-1} (M_i - M_N) \mathbf{J}_i\),
where \(\mathbf{v}\) is the coupled bulk velocity, \(M_i\) are molecular weights, and \(\rho\) is the
mixture density.

This material is designed to pair with flux materials such as
[`ADFluxFromGradientMaterial`](ADFluxFromGradientMaterial.md) that produce the individual species
fluxes.

## Example Input File

```text
[Materials]
  [molar_velocity]
    type = ADMolarAvgVelocityMaterial
    fluxes = 'J1 J2'
    molecular_weights = '2.0 16.0 28.0'
    density = rho
    molar_avg_velocity = u_m
    velocity = v
  []
[]
```

!syntax parameters /Materials/ADMolarAvgVelocityMaterial

!syntax inputs /Materials/ADMolarAvgVelocityMaterial

!syntax children /Materials/ADMolarAvgVelocityMaterial
