# Asteroid Mining Tech Demo: Trajectory Design Analysis

![MATLAB](https://img.shields.io/badge/MATLAB-R2024b-orange)
![Status](https://img.shields.io/badge/Status-Completed-green)
![Course](https://img.shields.io/badge/Course-Astrodynamics-blue)

## 📌 Abstract
This repository contains the solution for **Assessed Exercise 3** of the *Programming and Mathematics for Astrodynamics and Trajectory Design (2024-2025)* module.

The project involves a grid search for asteroid sample return opportunities for commercial nanosatellite demonstrators. The mission targets Near-Earth Asteroid (NEA) mining, specifically identifying **Asteroid 2014 WX202** as the optimal candidate. Transfer trajectories between early 2033 and the end of 2037 were analyzed to meet strict hyperbolic excess velocity and rendezvous maneuver constraints.

## 🎯 Mission Requirements
*   **Launch Window:** Jan 1, 2033 – Dec 31, 2037.
*   **Mission Profile:** Earth $\to$ Asteroid (Stay 2-6 months) $\to$ Earth.
*   **Constraints:**
    *   Departure/Arrival $v_{\infty} < 1.5$ km/s.
    *   Rendezvous maneuvers $\Delta v < 500$ m/s.

## ⚙️ Methodology

Due to the dataset size (~20,000 asteroids), a brute-force Lambert search was infeasible. A multi-step pruning process was implemented to filter candidates down to a manageable number before performing high-fidelity analysis.

### Process Overview
The selection strategy involved an initial pruning using the Shoemaker-Helin approximation, followed by a custom Figure of Merit (FoM) filter, and finally a rigorous Lambert arc scan.

![Methodology Block Diagram](assets/methodology_diagram.png)
*> Place your block diagram screenshot here (e.g., the flow chart from your PDF)*

---

### 🧮 Mathematical Proof: The Pruning Figure of Merit
To efficiently select the top 100 candidates, a custom **Figure of Merit (FoM)** was derived. This metric approximates the $\Delta v$ required to match the asteroid's orbit based on its eccentricity ($e$) and inclination ($i$), assuming low inclination/eccentricity and a circular Earth orbit.

#### I. Change Eccentricity ($\Delta v_1$)
The change in velocity required is the difference between the periapsis velocity ($v_p$) and the circular velocity ($v_c$).

$$
\Delta v_1 = |v_p - v_c|
$$

Using the *vis-viva* equation where $r = a$, the circular velocity is:
$$
v_c = \sqrt{\frac{\mu}{a}}
$$

The periapsis velocity at $r_p = a(1-e)$ is:
$$
v_p = \sqrt{\frac{\mu}{a}} \cdot \sqrt{\frac{1+e}{1-e}}
$$

Using a 1st order Taylor approximation for $e < 0.2$:
$$
\sqrt{\frac{1+e}{1-e}} \approx 1+e
$$

Substituting this back yields:
$$
\Delta v_1 \approx \left| \sqrt{\frac{\mu}{a}} \cdot (1+e) - \sqrt{\frac{\mu}{a}} \right| = v_c \cdot e
$$

#### II. Change of Inclination ($\Delta v_2$)
With no change of velocity magnitude (preserving semi-major axis $a$):

$$
\Delta v_2 = 2 v_c \cdot \sin(\Delta i / 2)
$$

#### III. Final Figure of Merit
Since the two maneuvers are applied perpendicularly ($\Delta v_1$ in-plane, $\Delta v_2 \perp$ plane), the total $\Delta v$ is:

$$
\Delta v_{tot} = \sqrt{(v_c \cdot e)^2 + \left( v_c \cdot [2 \cdot \sin(i/2)] \right)^2}
$$

Normalizing by $1/v_c$ gives the final FoM used in the code:

$$
\boxed{FoM = \sqrt{e^2 + [2 \cdot \sin(i/2)]^2}}
$$

---

## 🏆 Results

The analysis identified **Asteroid 2014 WX202** as the prime candidate. The recommended mission profile is as follows:

| Description | Date (MJD) | Date (Gregorian) | $\Delta v$ (km/s) |
| :--- | :--- | :--- | :--- |
| **Launch** | 12363.5 | Nov 2033 | - |
| **Arrival** | 12673.5 | Sep 2034 | 0.54 |
| **Departure** | 12733.5 | Nov 2034 | 0.54 |
| **Return** | 12923.5 | May 2035 | - |
| **Total ToF** | - | - | 560 Days |

### Porkchop Plots
Below are the porkchop plots generated for the four distinct burns of the mission.

<div align="center">
  <img src="assets/porkchop_dv1.png" width="45%" alt="Porkchop Plot dv1" />
  <img src="assets/porkchop_dv2.png" width="45%" alt="Porkchop Plot dv2" />
</div>
<div align="center">
  <img src="assets/porkchop_dv3.png" width="45%" alt="Porkchop Plot dv3" />
  <img src="assets/porkchop_dv4.png" width="45%" alt="Porkchop Plot dv4" />
</div>

*> Replace the above paths with the actual screenshots of your plots.*

## 📂 Repository Structure

```text
├── Code/
│   ├── main_script.m       # Main execution file
│   ├── lambert_solver.m    # Implementation of Izzo's algorithm
│   └── tools/              # Helper functions
├── Reports/
│   ├── Claudio_Ferrara_report.pdf  # Final Executive Report
│   └── fom.md              # Original Math Proof Markdown
└── assets/                 # Images for README
