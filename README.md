# Orbital Dynamics of Earth and Jupiter

This repository presents a numerical study of the orbital motion of Earth and Jupiter around the Sun. The project was originally completed as a final exam for a physics class in 2024 and has been cleaned and documented for professional presentation.

## 🌌 Overview

The goal of this project is to model and visualize the trajectories of Earth and Jupiter using classical Newtonian gravitation. The simulation integrates the equations of motion numerically and compares the resulting orbits with expected astronomical parameters.

The project includes:
- A two‑body and three‑body gravitational model  
- Numerical integration of orbital equations in C++  
- Visualization of orbital trajectories using Python  
- Discussion of orbital stability and perturbations  

## 🧠 Physics Background

The motion of each planet is governed by Newton’s law of gravitation:

$$
\vec{F} = -G \frac{m_1 m_2}{r^2} \hat{r}
$$

which leads to a system of second‑order differential equations for position and velocity.

The simulation uses:
- Gravitational constant $G$  
- Realistic planetary masses  
- Simplified initial positions and velocities  

## 🛠️ Technologies Used

- **C++** for numerical integration  
- **Python / Jupyter Notebook** for plotting (NumPy, Matplotlib)  
- **GitHub** for version control and documentation  

## 📁 Repository Structure

```text
src/
    orbit_simulation.cpp    # Numerical integration of the orbits
    plotting.ipynb          # Jupyter Notebook for visualization
    figures/                # Generated plots
README.md
LICENSE
```

## 📊 Example Output

The notebook generates orbital plots showing:
- Elliptical trajectories of Earth and Jupiter  
- Relative scales of the two orbits  
- Comparison with expected orbital periods  

## 🚀 How to Run

### 1. Compile and run the C++ simulation
```bash
g++ src/orbit_simulation.cpp -o orbit_simulation
./orbit_simulation
```

This produces the numerical data used for plotting.

### 2. Open the Jupyter Notebook
```bash
jupyter notebook src/plotting.ipynb
```

Running the notebook will generate the figures in `src/figures/`.

## 📚 About This Project

This work was originally completed as a final exam in a physics course in 2024. It has been updated and documented for clarity and professional presentation.
