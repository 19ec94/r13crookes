# r13 crooks

## Description

`r13 crooks` is a project that simulates the behavior of a **Crooks radiometer** using the **R13 moment equations**. This project contains source code that models the radiometer's performance across various geometries, including circle, rectangle, and diamond-shaped vane simulations. The radiometer's operation is governed by thermodynamic principles, and the R13 equations provide a framework for understanding its dynamics.

## Required Software

To run and interact with this project, you will need to install the following software:

- **Gmsh**: A 3D finite element mesh generator for scientific computing.
- **Docker**: A platform for developing, shipping, and running applications inside containers.
- **Apptainer**: A tool for running containerized applications with high performance and flexibility.
- **Python**: The programming language used for post-processing and calculation tasks.

## Project Structure

The repository contains the following folders and their respective content:

- **`fenicsR13`**: This folder contains the source code implementing the R13 moment equations and the simulation logic for the Crooks radiometer. This is adopted from Lambert's work.
- **`3d_crooks`**: This folder contains necessary input files and post-processing calculations. It includes simulations for various geometries, such as:
  - Circle-shaped vane
  - Rectangle-shaped vane
  - Diamond-shaped vane

## Getting Started

Follow the steps below to set up the project on your local machine.

### Prerequisites

Make sure that the required software listed above is installed.

### 1. Clone the Repository

Clone the repository to your local machine using the following command:

```bash
git clone https://github.com/19ec94/r13crooks.git
cd r13crooks
```

### 2. Set up the Simulation Environment

you can use Apptainer to run the containers. Ensure you have Apptainer installed and then run:

```bash
apptainer run fenicsr13.sif
```
