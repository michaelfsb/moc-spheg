# moc-spheg [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=michaelfsb/moc-spheg)

## Overview
This repository contains the code related to the research article titled "Control óptimo multifase de sistema fotovoltaico-hidrógeno asistido por red," which will be presented at the XX Reunión de Trabajo en Procesamiento de la Información y Control (RPIC 2023) from November 1 to November 3, 2023.

## Authors
- Michael Barbosa (Universidade Federal de Santa Catarina)
- José García Clúa (Instituto LEICI, Facultad de Ingeniería, UNLP)
- Gustavo de Andrade (Universidade Federal de Santa Catarina)
- Julio Normey-Rico (Universidade Federal de Santa Catarina)

## Abstract
Water electrolysis is being increasingly considered worldwide because of its potential to produce green hydrogen using renewable energy sources. This work focuses on the control of a PEM electrolyzer powered by photovoltaic solar energy, assisted by an electrical grid. To minimize this assistance, the application of numerical optimal control techniques is proposed, particularly the Hermite–Simpson direct collocation method. To meet the variable operating specifications of the electrolyzer, both single-phase and multiphase techniques have been suggested. These techniques predefine the temporal subintervals for switching dynamics and restrict the optimization variables for various conditions. Both frameworks are transformed into an nonlinear program (NLP) problem, which is solved using the IPOPT optimizer. The results obtained enable the analysis and comparison of improvements introduced not only in terms of NLP resolution time but also in the volume of hydrogen produced, utilization of renewable energy, and minimization of grid energy consumption.

## Dependencies
This code was developed and tested using the following software and versions:
- [MATLAB](https://www.mathworks.com/products/matlab.html) R2023a
- [CasADi](https://web.casadi.org/) v3.6.3

## How to Run
You can execute the code in two ways: locally or in MATLAB Online. In both cases, the first step is to install CasADi.

### Local Execution:
1. Install CasADi by following the instructions provided on the [CasADi website](https://web.casadi.org/get/).
2. Obtain the code by either:
   - Cloning this repository using the command:
     ```
     git clone https://github.com/michaelfsb/moc-spheg.git
     ```

### MATLAB Online Execution:
1. Install CasADi by following the instructions provided on the [CasADi website](https://web.casadi.org/get/).
2. Use this link to clone this repository to your MATLAB Online environment:<br> 
     https://matlab.mathworks.com/open/github/v1?repo=michaelfsb/moc-spheg


Afterward, you can access the "src" folder and run the following scripts: "single_phase.m", "multi_phase_1.m" or "multi_phase_2.m" to obtain the results presented in the article.
