# FluidSim
A simple real-time 2D fluid simulation using Smoothed Particle Hydrodynamics (SPH) with OpenGL compute shaders.

![FluidSim Demo](assets/fluid_sim_demo.gif)

## Overview
This project implements an SPH particle simulator in C++ using GLFW for windowing, GLAD for OpenGL function loading, and GLSL compute shaders for physics updates. Particles are rendered as instanced quads, and spatial hashing is used to group particles into grid cells for neighbor queries.

## Features
- GPU-accelerated physics using OpenGL compute shaders
- Density, pressure, viscosity, and gravity calculations
- Particle collision handling with screen boundaries
- Mouse interaction with left/right button forces
- Real-time rendering using instanced geometry

## Build
The project is configured to compile with `g++` on Windows.

Example build command:
```sh
cd c:\Users\lawpe\Dev\FluidSim
g++ -g -Iinclude -Llib src\*.cpp src\glad.c include\imgui\*.cpp -lglfw3 -lopengl32 -lgdi32 -o fluidsim.exe
```

## Run
After building, run the executable:
```sh
.\fluidsim.exe
```

## Controls
- `Space`: pause/resume simulation
- `.` (period): step one frame when paused
- `,` (comma): step up to 10 frames when paused
- `R`: reset particles to the initial grid
- Left mouse button: interact with particles
- Right mouse button: alternate mouse interaction mode

## Project Structure
- `src/main.cpp` - application entry point, OpenGL setup, input callbacks, compute shader dispatch, and rendering loop
- `src/particle.hpp` - particle data structure, collision logic, and utility methods
- `include/shaders/` - GLSL shader sources for rendering and compute-based fluid simulation
- `include/` - third-party headers for GLAD, GLFW, GLM, and ImGui

## Notes
The simulation expects OpenGL 4.6 compatibility for compute shader support and uses Shader Storage Buffer Objects (SSBOs) to transfer particle data between CPU and GPU.
