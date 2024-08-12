# HFQsim
A Simulator of Earthquakes and aseismic slip on a Heterogeneous strike-slip Fault (HFQsim) with static/kinetic friction and temperature-dependent creep
## Abstract

A Simulator of Earthquakes and aseismic slip on a Heterogeneous strike-slip Fault (HFQsim) with static/kinetic friction and temperature-dependent creep

Xiaoyu Zhou1, Yehuda Ben-Zion1,2 
1. University of Southern California, Department of Earth Science 
2. Statewide California Earthquake Center


We develop an earthquake simulator to study the partitioning of seismic/aseismic slip and dynamics of Earthquakes on a Heterogeneous strike-slip Fault (HFQsim) using a generalized model of a discrete fault governed by static/dynamic friction and creep in an elastic half-space. Previous versions of the simulator were shown to produce various realistic seismicity patterns (e.g., frequency-magnitude event statistics, hypocenter and slip distributions, temporal occurrence) using friction levels and creep properties that vary in space but are fixed in time. The new simulator incorporates frictional heat generation by earthquake slip leading to temperature rises, subsequent diffusion cooling into the half space, and time-dependent creep on the fault. The model assumes a power law dependence of creep velocity on the local shear stress, with temperature-dependent coefficients based on the Arrhenius equation. Temperature rises due to seismic slip produce increased aseismic slip, which can lead to further stress concentrations, aftershocks, and heat generation in a feedback loop. The partitioning of seismic/aseismic slip and space-time evolution of seismicity are strongly affected by the temperature changes on the fault. The results are also affected significantly by the difference between the static and kinetic friction levels. The model produces realistic spatio-temporal distribution of seismicity, transient aseismic slip patterns, foreshock-mainshock-aftershock sequences, and a bimodal distribution of earthquakes with background and clustered events similar to observations. The HFQsim may be used to clarify relations between fault properties and different features of seismicity and aseismic slip, and to improve the understanding of failure patterns preceding large earthquakes. The HFQsim is an open-source code in the process of being hosted by the science gateway Quakeworx. 
