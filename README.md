#MomentumScale 

This package is able to decompose the momentum scale biases in a global scale $`\varepsilon_{s}`$ and a residual Z $`\varepsilon_{dz}`$ (or alternatively radial) scale based on the following equations:

```math
{m'}_{\mu\mu}^2/m_{\mu\mu}^2 - 1\approx + 2 A^{+}_s \varepsilon_{s}(\eta^+,\phi^+) + 2 A^{-}_s \varepsilon_{s}(\eta^-,\phi^-) + 2 A^{+}_z \varepsilon_{dz}(\eta^+,\phi^+) 	+ 2 A^{-}_z \varepsilon_{dz}(\eta^-,\phi^-) \\ 
```

where $`A^{\pm}_s =  E^{\pm}E^{\mp}\left( \pmb{\beta}^{\pm} - \pmb{\beta}^{\mp} \right )^2 /m_{\mu\mu}^2 `$
and  $`A^{\pm}_z =  E^{\pm}E^{\mp}\left[ \left(\beta_{\mathrm{z}}^{\pm} \right)^2 - \pmb{\beta}_{\mathrm{z}}^{\mp} \cdot \pmb{\beta}_{\mathrm{z}}^{\pm}\right] /m_{\mu\mu}^2  `$.

And the parameters are constrained such that:
$` E^+E^-\left( \pmb{\beta}^+ - \pmb{\beta}^- \right)^2 =  m^2_{\mu\mu}  - m^2_\mu \left(E^+ + E^- \right)^2/E^+E^-  `$  
and that 
$` E^+E^-\left( \left( \beta^+\right)^2 - \pmb{\beta}^+ \cdot \pmb{\beta}^- \right) =  m^2_{\mu\mu}/2 - m_{\mu}^2 (1+E^{-}/E^{+})  `$
