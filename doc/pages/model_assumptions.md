---
layout: default
title: Model assumptions
mathjax: true
---

# Model assumptions

## Model geometry

QDYN handles the following geometries:
- **0D fault**, a spring-block model
- **Finite 1D fault** embedded in a 2D unbounded elastic medium. The fault is straight and actually infinitely long, but its frictional part is confined to a segment of finite length $L$, $x \in \left[-L/2, L/2 \right]$.
- **Periodic 1D fault** embedded in a 2D unbounded elastic medium. The fault is infinitely long but has a spatially periodic distribution of friction parameters, normal stress and slip with spatial period $L$. The modeled segment is $x \in \left[-L/2, L/2 \right]$. Other variants: 
  - Fault bisecting an elastic slab of uniform thickness, with uniform steady displacement applied along the boundaries of the slab. 
  - Fault bisecting a low rigidity layer of uniform thickness embedded in a stiffer elastic medium.
- Vertical, anti-plane 1D fault in a half space.
- 1.5D fault model in which we account approximately for a characteristic fault width W, representing for instance the seismogenic width of a 2D fault. Slip is approximated by a sinusoidal pattern over a length W in the out-of-plane dimension. The computational cost of this reduced-dimensionality problem is the same as that of 1D-fault simulations.
- 2D fault embedded in a 3D elastic space or half-space. The fault surface has fixed strike, but possibly depth-dependent dip. The fault is infinite but only a finite area is frictional.



## Boundary conditions

Spring-block models are loaded by an imposed load-point velocity. On continuum faults, the frictional segment is loaded by slip imposed along the remaining, non-frictional part of the fault. In all cases, the imposed loading is composed of a steady velocity and an oscillatory component. 

In the quasi-dynamic approximation adopted in QDYN, fault stresses are the sum of static elastic stresses induced by slip and a radiation damping stress. The latter approximates the effect of wave radiation: it represents exactly the stresses induced by waves radiated in the direction normal to the fault but not the
complete elastodynamic stresses.



## Fault rheology

The behavior of the model fault in response to the imposed boundary conditions is governed by the fault rheology. QDYN offers two classes of fault rheology: *rate-and-state friction* and the *Chen-Niemeijer-Spiers* microphysical model.

### Rate-and-state friction laws

The friction coefficient is governed by one of the following rate-and-state friction laws:

- Conventional rate-and-state: $\mu (V, \theta) = \mu^* + a \ln \left( \frac{V}{V^\*} \right) + b \ln \left( \frac{\theta V^\*}{D_c} \right)$
- Regularized rate-and-state friction: $\mu (V, \theta) = a \sinh^{-1} \left( \frac{V}{2V^*} \right) \exp \left( \frac{1}{a} \left[\mu^\* + b \log \left( \frac{\theta V^\*}{D_c} \right) \right] \right)$
- Rate-and-state with cut-offs: $\mu (V, \theta) = \mu^* - a \ln \left( 1 + \frac{V_1}{V} \right) + b \ln \left( 1+ \frac{\theta V_2}{D_c} \right)$

The state variable follows one of the following evolution equations:

- Aging law: $\frac{\mathrm{d} \theta}{\mathrm{d} t} = 1 - \frac{\theta V}{D_c}$
- Slip law: $\frac{\mathrm{d} \theta}{\mathrm{d} t} = - \frac{\theta V}{D_c} \ln \left( \frac{\theta V}{D_c} \right)$

All frictional parameters can be spatially heterogeneous.



### Microphysically-based friction laws

As an alternative to rate-and-state friction, simulations can be run using a microphysical formulation for the fault rheology. The fault mechanics model is based on the Chen-Niemeijer-Spiers (CNS) model (*Niemeijer & Spiers*, [2007](https://doi.org/10.1029/2007JB005008); *Chen & Spiers*, [2016](https://doi.org/10.1002/2016JB013470)). The implementation into QDYN is detailed in *Van den Ende et al.* ([2018](https://doi.org/10.1016/j.tecto.2017.11.040)). Conceptually, the CNS model is based on the interplay between dilatant granular flow and non-dilatant ductile creep, from which both velocity-strengthening and velocity-weakening behaviour may emerge.