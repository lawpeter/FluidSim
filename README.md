# FluidSim

An interactive GPU-accelerated 2D fluid simulation using Smoothed Particle Hydrodynamics (SPH),
built in C++ with OpenGL compute shaders for GPU-accelerated physics.

## Overview

This project uses a hybrid CPU/GPU pipeline: particle state is read back to
the CPU for spatial-grid construction, index data is uploaded, and SPH work is
dispatched through GLSL compute shaders using OpenGL Shader Storage Buffer Objects (SSBOs). Particles
are color-coded by velocity, giving a vivid visual sense of fluid behavior —
from tight, pressurized splashes to slow gas-like expansion.

Built collaboratively as both a fluid simulation and a learning exercise in
GPGPU programming. We used Sebastian Lague's fluid simulation video as
inspiration, read the same papers he referenced, and derived the SPH math
ourselves. His implementation uses HLSL; we chose GLSL deliberately to learn
a different shader pipeline.

**Reference:** [Sebastian Lague — Coding Adventure: Fluid Simulation](https://youtu.be/rSKMYc1CQHE?si=b2XPkc6CLPTXR4fl)

## Features

- GPU-accelerated SPH physics via GLSL compute shaders
- A uniform spatial grid restricts neighbor searches to nearby cells; work depends on cell occupancy
- Density, pressure, viscosity, and gravity force calculations
- Velocity-based particle color coding
- Mouse interaction — push or pull particles with configurable strength and radius
- Adjustable simulation parameters: gravity, target density, mouse force
- Particle collision response at screen boundaries
- Interactive rendering using instanced quad geometry

## Requirements

- OpenGL 4.6 (compute shader and SSBO support required)
- Windows (main branch) or Apple Silicon macOS (mac branch)

## Build

### Windows (main branch)

From the repo root:

```bash
g++ -g -Iinclude -Llib src\*.cpp src\glad.c include\imgui\*.cpp \
    -lglfw3 -lopengl32 -lgdi32 -o fluidsim.exe
```

### macOS — Apple Silicon (mac branch)

Switch to the mac branch, then from the repo root:

```bash
/usr/bin/g++ -g -Iinclude -Llib src/*.cpp src/glad.c include/imgui/*.cpp \
    -lglfw3 \
    -framework OpenGL \
    -framework Cocoa \
    -framework IOKit \
    -framework CoreVideo \
    -framework CoreFoundation \
    -o fluidsim.o
```

## Run

**Windows:**
```bash
.\fluidsim.exe
```

**macOS:**
```bash
./fluidsim.o
```

## Controls

| Input | Action |
|---|---|
| `Space` | Pause / resume |
| `.` | Step one frame (while paused) |
| `,` | Step up to 10 frames (while paused) |
| `R` | Reset particles to initial grid |
| Left mouse button | Push particles away |
| Right mouse button | Pull particles toward cursor |

Simulation parameters (gravity, target density, mouse strength and radius)
are adjustable via the ImGui panel at runtime.

## Project Structure
src/
main.cpp          — OpenGL setup, input handling, compute dispatch, render loop
particle.hpp      — Particle data structure, collision logic, utility methods
include/
shaders/          — GLSL shader sources (rendering + compute)
imgui/            — ImGui headers
glad.h / glfw /   — OpenGL loader and windowing
glm/              — Math library

## Status and Performance

Complete; no active development is planned. Peter observed approximately 290,000 particles at around 35 FPS in repeated runs on an NVIDIA RTX 5070 Ti / AMD Ryzen 9 9950X desktop. This was an informal observation, not a controlled benchmark. Frame rate does not establish synchronization between simulation time and wall-clock time.

## Limitations

- Behavior was evaluated qualitatively; no quantitative comparison against analytical, experimental, or validated CFD reference cases was performed.
- CPU grid construction and readback/upload introduce CPU/GPU synchronization overhead.

- Simulation becomes unstable and may crash above ~300,000 particles;
  root cause not yet identified
- No 3D support; extending to 3D SPH is a natural next step
- Further GPU-side optimizations (e.g. prefix sum compaction, better grid
  boundary handling) were scoped out but not implemented

## Credits

Developed collaboratively by [lawpeter](https://github.com/lawpeter) and
[owendpoole](https://github.com/owendpoole).

SPH math derived from the following papers, as referenced by Sebastian Lague:
- Müller et al., *Particle-Based Fluid Simulation for Interactive Applications* (2003)
- Monaghan, *Smoothed Particle Hydrodynamics* (1992)
