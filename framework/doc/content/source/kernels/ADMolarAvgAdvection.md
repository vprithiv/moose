# ADMolarAvgAdvection

!syntax description /Kernels/ADMolarAvgAdvection

## Overview

`ADMolarAvgAdvection` applies advection using a molar-average velocity provided by a material
property. The kernel evaluates the contribution
\(R_i = \left( \mathbf{u}_\mathrm{m} \cdot \nabla u_h, \psi_i \right)\) for the advected species
variable, where \(\mathbf{u}_\mathrm{m}\) is typically supplied by
[`ADMolarAvgVelocityMaterial`](../materials/ADMolarAvgVelocityMaterial.md).

The kernel requires the name of the molar-average velocity material property and multiplies that
vector with the gradient of the primary variable at each quadrature point.

## Example Input File

```text
[Kernels]
  [molar_advect]
    type = ADMolarAvgAdvection
    variable = c1
    molar_avg_velocity = u_m
  []
[]
```

!syntax parameters /Kernels/ADMolarAvgAdvection

!syntax inputs /Kernels/ADMolarAvgAdvection

!syntax children /Kernels/ADMolarAvgAdvection
