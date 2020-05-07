# dmrtml
The Dense Media Radiative Transfer - Multi Layers model

## What is DMRT-ML ?

DMRT-ML is a physical model used to compute the thermal microwave emission of a given snowpack for passive microwave remote sensing applications. The model is based on the Dense Media Radiative Transfer Theory (Tsang, 1992 and Tsang and Kong 2001) and accurately solves the radiative transfer equation using the Discrete Ordinate Method (DISORT, Jin, 1994).

The snowpack is modeled as a stack of horizontal layers of snow and an optional underlying interface representing the soil or the ice. The atmospheric downwelling contribution can be taken into account.

DMRT-ML is designed to work for most snow-covered surfaces, and can model dry or wet snowpacks over soil (e.g. Alpine or Arctic seasonal snow), over ice (e.g. on ice-sheet or lake) and soon over sea-ice.

DMRT-ML has been developed at LGGE (now IGE) and is freely distributed under an open source license. Contributions from the remote sensing and snow community would be appreciated to further evaluate, enhance and document DMRT-ML.

The model is written in Fortran 90 and (optionally) interfaced to Python (version 2.5 and higher). It is thus fast and easy to use. The code is known to compile and run on Linux and Windows and should run smoothly on other systems.

Detailed description and analysis of several simulations are given in the reference paper

## How to run DMRT-ML ?

See the examples in the example directoty.

## Main and contributing authors :

Ghislain Picard1, Ludovic Brucker1*, Alexandre Roy2, Florent Dupont1,2

- 1 Insitut des Geosciences de l'Environnement (IGE)
54 rue Molière - Domaine Universitaire - BP 96
38402 St Martin d'Hères Cedex, FRANCE
- 2 Centre d'applications et de recherches en télédétection (CARTEL)
Université de Sherbrooke
2500 Bd Université
Sherbrooke, Qc J1K 2R1 CANADA
- (*) now at: NASA GSFC, Cryospheric Sciences Lab., code 615 Greenbelt, MD 20771 U.S.A. 

