# The problem



Two adherents with equal length X and Y, defined by the user in ["simulation"]["lenX"] and ["simulation"]["lenY"] respectivelly, being welded in through-the-thickness direction.

Their initial temperature is defined by the user in ["Materials"]["Material1"]["Domain Initial Temperature"] and ["Materials"]["Material2"]["Domain Initial Temperature"] respectivelly.

In case of isothermal simulation, i.e forcing thier interface to be at constant temperature, the user defines the interface temperature in ["Processing Parameters"]["Temperature"].

In case of anisothermal simulation with constant power density, the used define the value (W/m^3) in ["Processing Parameters"]["Temperature"].

At this point, thermal conductivity (W/(mºC)), Specific Heat (J/(J/KgºC)) and density (Kg/m^3) do not change with temperature nor time. The user define those values in ["Materials"].

The user also define the number of elements in X and Y in ["simulation"]["Number of Elements X"] and ["simulation"]["Number of Elements Y"] respectivelly.

The number of steps and step time are defined in ["simulation"]["Step Time"] and ["simulation"]["Number Time Steps"] respectivelly. Notice that large step time may not converge and an error message will popup.

Finally, the code generates Dic and T plots at every interval, defined in ["Plot"]["plot interval"].



```yaml
Problem Type:
    Type: TwoPlates
Processing Parameters:
    Temperature: 648
    Pressure: 0.1
    Power Density: 2613203000
Materials:
  Material1:
    Thermal Conductivity X: 0.65
    Thermal Conductivity Y: 0.65
    Density: 1598
    Cp: 930
    Domain Initial Temperature: 293.
  Material2:
    Thermal Conductivity X: 0.65
    Thermal Conductivity Y: 0.65
    Density: 1598
    Cp: 930
    Domain Initial Temperature: 293.

Simulation:
  Time Step: 0.0001
  Number Time Steps: 10001
  lenX: 0.2
  lenY: 0.005248
  Number of Elements X: 200
  Number of Elements Y: 128

Plot:
  folder1: "./output/Temperature/"
  folder2: "./output/Dic/"
  figure temperature name: Temperature
  figure dic name: Dic
  Color Interpolation: 50
  Color Map: "inferno"
  plot interval: 500
  
Animation:
    name: temperature
```


# Getting Started

### Mesh

A class named "Mesh Two Plates" is defined to create the mesh and assign Initial Temperature, Thermal Diffusivity in X and Thermal Diffusivity in Y in each element.
### Model
A class named "Heat Transfer" solves the heat transfer problem by central differences.
A class named "Intimate Contact" calculates the viscosity in function of the temperature and the evolution of the intimate contact.
### Plot and Animation
A class named "PlotsTwoPlates" is defined to generate figures of the thermal history in spaced intervals. In the end, a .GIF is created from those figures.
### Solver
A class named "SolvesTwoPlates" is defined to create a forward step analysis, solving the heat transfer model and plotting the figures when convenient.