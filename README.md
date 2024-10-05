# Polydisperse molecular dynamics in two dimensions

This code is used to simulate a polydisperse mixture of hard-like particles
in two dimensions. It uses molecular dynamics to propagate the positions
of the particles.

## Features

A short list of features:

- Uses cell lists from `CellListMap.jl`; the parallelization is turned off by default.
- The integrator is velocity Verlet.
- The polydisperse mixture is obtained by sampling values from a normal distribution and assigning these values to the diameter of the particles. There is a safeguard to avoid having diameters larger than 1.
- The temperature is controlled using the Bussi-Donadio-Parrinello thermostat.
