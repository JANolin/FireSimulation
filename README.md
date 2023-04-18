
# Simplified 2D Physically Based Fire Simulation

## Description
This is an implementation of a simple fire animation realised for COMP 559 extended off the code for assignment 3. This is largely based off Real-Time Fluid Dynamics for Games by Jos Stam and Realistic 2D Fire in Real-Time by Odd Erik Gundersen et al.

## Installation
This code requires the same dependencies as A3, which I have provided in the repo under thirdparty. Here are the instruction from Prof. Kry on how to set them up. https://github.com/paulkry/comp559cppSetup


## Usage
There is two menus, one for the control of parameters and the other for visualization. In the parameter menu, you can switch between temperature source and fuel source to add a source to the corresponding density field. In the visualization menu, there is three different rendering possible: temperature which only takes into account the temperature (useful to start a fire), fuel which is a quick hack to check the amount of fuel in the system and fire which is the visualization of fire as described in the paper. All visualization use the black-body radiation model, but only the fire visualization is accurate for the rendering: the other two are just for setting up the fire reaction.

## Simulating Fire
Creating a self-sustaining chain-reaction is incredibly hard. However, if you keep an eye on the dissipation parameter and adjust it as the reaction progresses, you can approximate the effect to prevent the system from exploding, or the flame from dying out. To start a flame, I first add a couple of heat source at max srcTemp. When I reach a temperature over the reactionTemperature, I add a fuel source at max srcTemp (here, srcTemp refers to the source value). If a good balance in the parameters is achieved, a flame should be achieved. Don't forget to switch to fire visualization after the fuel source is added.

## Parameters
Here is a list of all the new controllable parameters and their effect:

* dissipation: controls the amount of dissipation in each corresponding density field. Useful to prevent explosion.

* gravity: straightforward

* reactionTemperature: controls the temperature in each grid square required to determine if combustion is happening in that square

* expansion: scales the horizontal expansion you would normally get from a flame. Approximates the core of the flame to be at the horizontal center of the grid

* pctFuelBurn: Percentage of fuel burned by combustion in each square after one second.

* stochiometricMixture: stochiometric amount of oxygen needed for one unit of fuel to burn. A value of 2 would halves the reaction since you would need twice the amound of oxygen in each grid; hence only half of the available fuel would burn.

* combustionHeat: amount of heat created by the combustion reaction. Careful, as this can easily make your system explode (too much heat).

* vorticityStrength: scales the vorticity, which adds the turbulent nature to the flame.

And in the visualization menu, we have:

* mapIntensity: Scales the intensity of the light emitted by the reaction. I would recommend to lower it when you try to get the fire started to get a better view of the flame. High values = more obscure reaction