# Vicsek Model Simulation in MATLAB

## Overview
This code simulates a basic Vicsek model in which the collective motion of particles is influenced by alignment with their neighbors. In addition, it analyzes the order parameter (average normalized velocity) at different levels of noise. 

## Usage

### 1. Initialization of Parameters
This section sets the simulation parameters and initiaizes random particle positions and angles of movement. 

### 2. Main Simulation Loop
This section simulates particle movements in a periodic boundary environment. At each time step, it (a) updates each particle's angle of movement, (b) updates each particle's position, and (c) plots the data. In addition, this section is used to calculate the order parameter and to check for system stability. 

### 3. Order Parameter Analysis
  #### i. As a function of time
  This section plots the order parameter over time for the last stored simulation. 

  #### ii. As a function of noise
  This section plots the steady-state order parameter for the last stored simulation as a function of noise level (for different combinations of particle number and box size).

## Results
- Scatter plot figure
- Optional Video: `Vicsek_Model_Simulation_condition_X.mp4`
  - Visualization of particle positions over time.
- Plot: `Va_plot_over_time.png`
  - Order parameter evolution over time.
- Plot: `Steady_State_Velocity_vs_Eta_Multiple.png`
  - Steady-state average velocity as a function of noise level for different parameter combinations.

