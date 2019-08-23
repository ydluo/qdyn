---
layout: default
title: Model assumptions
mathjax: true
---

# Model assumptions

## Model geometry

QDYN handles the following geometries:
- **0D fault**: a spring-block model

- **Finite 1D fault**: a straight fault (a 1D line) embedded in a 2D unbounded elastic medium. Slip is antiplane, i.e. in the direction normal to the 2D plane. The fault is infinitely long, but friction boundary conditions are applied only on a segment of finite length $L$, $x \in \left[-L/2, L/2 \right]$. Outside this segment, constant slip velocity is prescribed.

- **Periodic 1D fault**: same as above, except that here friction is applied all along the infinite fault, and the spatial distributions of friction parameters, normal stress and slip are assumed to be periodic with spatial period $L$. The following variants are also available: 
  - Fault bisecting an elastic slab of uniform half-thickness $H$, with uniform, steady displacement applied along the boundaries of the slab located at a distance $H$ from the fault.
  - Fault bisecting a low rigidity layer of uniform half-thickness $H$ embedded in a stiffer elastic medium.
  
- Vertical, anti-plane, straight 1D fault in a 2D half space.

- 1.5D fault model (2.5D elastic medium) that accounts approximately for a finite fault size $W$ in the out-of-plane direction (normal to the strike direction $x$), assuming slip either tapers off to zero at the $W$-edge or steadily creeps beyond $W$. This represents for instance a map-view model of a vertical fault with finite seismogenic depth $W$. For more details and a derivation, see appendix A.2 of *Luo and Ampuero* ([2017a](http://dx.doi.org/10.1016/j.tecto.2017.11.006)). The computational cost of this reduced-dimensionality model is the same as that of 1D-fault simulations.

- 2D fault embedded in a 3D elastic space or half-space. The fault surface has fixed strike, but possibly depth-dependent dip. The fault is infinite but friction is applied only over a finite area.

  

## Boundary conditions

Spring-block models are loaded by an imposed displacement applied at a load-point connected to the rigid block by an elastic spring. 

On non-periodic, infinite faults embedded in a continuum elastic medium, friction boundary conditions are applied over a fault segment of finite size. This frictional segment is loaded by slip imposed along the remaining part of the fault. 

On periodic faults embedded in an unbounded elastic medium, the fault is loaded by slip prescribed uniformly on the fault portions beyond the frictional width $W$.

On periodic faults embedded in a bounded elastic medium of half-thickness $H$, the fault is loaded by displacement prescribed uniformly along the boundaries of the medium located at a distance H from the fault.

In all cases, the imposed slip or displacement usually has constant velocity. An additional oscillatory component is also possible.

On the frictional portion of the fault, the shear strength equals the normal stress times the friction coefficient. Faults governed by rate-and-state friction are always slipping (although very slowly when locked during interseismic periods) and their shear stress is always equal to their frictional strength. 

In the quasi-dynamic approximation adopted in QDYN, fault stresses are the sum of static elastic stresses induced by slip and a radiation damping stress. The latter approximates the effect of wave radiation: it represents exactly the stresses induced by waves radiated in the direction normal to the fault, but not the complete elastodynamic stresses.



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



### Thermal pressurisation

In addition to the (low-velocity) constitutive relations given by the rate-and-state and CNS models, QDYN includes high-velocity dynamic weakening through thermal pressurisation. The thermal pressurisation implementation is based on the spectral method of *Noda & Lapusta* ([2010](https://doi.org/10.1029/2010JB007780)), which solves the coupled differential equations for the diffusion of fluids and heat in the spectral domain. Diffusion only occurs normal to the fault plane, and it is assumed that the transport properties of the medium are homogeneous in the fault-normal direction and constant in time. In the fault-parallel direction, the transport properties may be spatially heterogeneous.