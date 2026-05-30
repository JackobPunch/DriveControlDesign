# DC Drive Cascaded Control System

University project for the *Electrical Drives* course at AGH University of Science and Technology (Kraków), Electrical Engineering. Design and simulation of a cascaded speed-current control structure for a separately excited DC motor.

## Motor Parameters

| P_N [kW] | U_N [V] | I_N [A] | n_N [RPM] | R_t [Ω] | L_t [mH] | J_S [kg·m²] |
|----------|---------|---------|-----------|---------|---------|------------|
| 17 | 230 | 85 | 700 | 0.253 | 1.9 | 0.75 |

## Scope

- **Mathematical model** — state-space and transfer function representation of the separately excited DC motor with thyristor converter
- **Cascaded control structure** — inner current PI controller (modular criterion) + outer speed PI controller (symmetry criterion)
- **Controller tuning** — current PI: G_RI(s) = (0.0085s+1)/0.6296s · speed PI: G_Rω(s) = (1.568s+10.89)/0.144s
- **Simulink implementation** — full closed-loop model with anti-windup saturation blocks
- **Startup simulations** — three scenarios: no-load then impact load, nominal active torque, nominal passive torque
- **Stability analysis** — Bode and Nyquist diagrams for the closed-loop system

## Simulation Results

Step response without cascaded control shows armature current and its derivative exceeding safe limits — confirming the need for the cascaded structure. With both PI controllers active, all three startup scenarios keep current within the 153 A limit and speed stabilises at the nominal 73.3 rad/s with minimal overshoot.

![Simulink model](https://github.com/JackobPunch/DriveControlDesign/blob/main/zrzutSimulink.png)

## Tools

MATLAB Simulink

## Authors

Jakub Cios, Maciej Duda — Electrical Engineering, AGH University of Science and Technology, Kraków
