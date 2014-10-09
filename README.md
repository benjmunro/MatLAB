MATLAB
======

> Various MATLAB scripts and functions, mainly written for Mechanical Engineering classes at the University of Canterbury

### Computational Fluid Mechanics


##### FVM Simulation of a Potential Flow around a 90° Corner

<p align="center">
  <img src="https://github.com/benjmunro/MatLAB/blob/master/Computational%20Fluid%20Mechanics/Stream_lines.png" alt="Stream lines"/>
</p>


##### Finite Volume Method Examples

Finite volume method codes for diffusion and convection-diffusion problems, examples mainly from the book An Itroduction to Computational Fluid Dynamics: The Finite Volume Method by Versteeg and Malalasekera.


**V & M Example 4.1**
Discription : Solve for temperature along a rod wih 1D heat conduction, no internal source or heat loss and specified end temperatures.

*OneDim_Heat_Conduction.m* -- solves Problem 4.1

**V & M Example 4.2**
Discription : Solve for temperature along a rod wih 1D heat conduction, with source, no heat loss and specified end temperatures.

*OneD_Heat_Conduction_in_a_Rod_with_Source.m* -- solves Problem 5.2

**V & M Example 4.3**
Discription : Cooling of a circular fin by means of heat convection along its length

*One_D_Rod_Cooling.m* -- Solves problem 4.3

**V & M Example 5.1**
Discription : A property is transported by means of convection and diffusion through a one dimensional domain

*Convetion_and_Diffusion_CD.m*  -- Solves Problem 5.1 with the Central Differencing Scheme
*Convetion_and_Diffusion_HYB.m* -- Solves Problem 5.1 with the Hybrid Differencing Scheme
*Convetion_and_Diffusion_UP.m* -- Solves Problem 5.1 with the Upwind Differencing Scheme

**V & M Example 8.1**
Discription : Solve for temperature along a thin plate wih tranisent 1D heat conduction,no internal source or heat loss and specified end temperatures.

*Unsteady_Heat_Transfer.m* -- Solves problem 8.1 with fully implict, fully explict or Crank Nickelson method 

### Vibrations and Dynamics
Methods for analyzing vibrations in mechanical systems.

Discrete Systems (Analysis)
  1. Time Domain Methods – Numerical Integration Methods
  2. Frequency Domain Methods
  3. Modal Domain Methods

→All methods and models will be developed for SDOF systems and applied to MDOF systems. Data for analysis is included, from the 2010 Canterbury Earth Quake. 

<p align="center">
  <img src="https://github.com/benjmunro/MatLAB/blob/master/Vibrations/Modal%20Displacement.jpg" alt="Stream lines"/>
</p>


##### Time Domain
*Euler_script.m* -- SINGLE degree of freedom system
*Newmark_Beta_MDOF_Script.m* -- mutiple degrees of fredom
*Newmark_Beta_SDOF_Script_ts.m* -- Single degree of fredom

##### Freqency Domain
*Freqency Domain Analysis-MDOF.m*  -- mutiple degrees of fredom
*Freqency_Domain_Analysis_SDOF.m*  -- Single degree of fredom

These scripts takes an input force, transfers it to the freqency domian using fft. The script calculates the trasfer function which is then mutiplied by the fft input force to get the response in the freqency domian which is then transfered back to the time domain using ifft

 *Modal_Domain_Analysis.m*
 
 This scripts takes an input force, transfers it to the modal domian. Returns Mode Shapes, Mode response and time response.  
 
 *DFT.m*
 
 Discrete Fourier Transform. Returns the real and imaginary parts of the F in the frequency domain
 
 *IDFT.m*
 
 Inverse Discrete Fourier Transform.Returns the time domain form from the real and imaginary parts of the F 
 in the frequency domain

*Euler.m*

Performs Euler Intergration, SDOF
 
*Newmark_B_SDOF.m* 
*Newmark_B_MDOF.m*

Newmrk Beta integration. Integrates a MDOF/SDOF system with mass Matrix "M", stiffness Matrix "K" and damping coeffiecient Matrix "D", when subjected to an external load f(t). Returns the displacement, velocity and acceleration of the system with respect to an inertial frame of reference.

### Controls and Parameter ID


<p align="center">
  <img src="https://github.com/benjmunro/MatLAB/blob/master/Controls/MC%20alalysis.png" alt="Stream lines"/>
</p>


### Computational Mathematics
