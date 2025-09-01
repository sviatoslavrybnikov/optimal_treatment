# Optimal Treatment Strategy
This repository contains the Octave/MATLAB code used in the study:

**"Integrating gene drives with established pest controls to mitigate spillover risk"**  
by *Sviatoslav R. Rybnikov, Adam Lampert, and Gili Greenbaum*.

---

## Overview
The code identifies the **optimal treatment strategy** involving:
- gene drive deployment  
- pesticide application  
- sterile male release  

The code employs dynamic programming to minimize the **total cost of eradication**, including:
- costs due to spillover risk
- expenditures on pesticide application and/or sterile male release

The code requires, as **input parameters**:
- gene drive configuration (conversion rate, fitness cost, and dominance)
- per-unit costs of possible interventions (pesticides and sterile males) relative to spillover risk
- initial conditions (population density and gene drive frequency)

---

## Requirements
- Octave 7.1 or later
