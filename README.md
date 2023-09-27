# ViAGRA
Tidal response of a Viscoelastic, Auto-GRAvitating body algorithm

This is a preliminary (but working!) version of VIAGRA (VIsco-elAstic response for an auto-GRAvitating body).

The code uses the equation of visco-elastic deformation for an auto-gravitating body under tides from a primary body.
It assumes that the studied body (typically the secondary) is locked and is restricted to degree 2 tides for now.
Models can have solid and liquid layers, each having their own set of equations being solved using matrix propagation.

To use it:
- Download all files and folders.
- First, add to the root a .txt file, tab separated, with the interior model to be used.
Columns must be radius (m), density (kg.m-3), first Lame parameter lambda (Pa), second Lame parameter mu (Pa) (aka shear modulus), dynamic viscosity (Pa.s-1)
- Add then a .m file containing the orbital informations of the system.
Parameters must have the same names as below:
M1: mass of the primary (kg)
a: semimajor axis of the studied body orbit around its primary (m)
e: orbital eccentricity
per: orbital period, used to derive the forcing frequency of the system. Can be changed to study Love numbers frequency dependency.

Obliquity tides are not yet implemented. When they are, obliquity will be required in this file.

- Fill the required parts in the Main.m file according to the previous steps and the wanted simulation.
Mostly as wanted in the EDITABLE section. MAIN CALL section typically should not be modified; proceed with caution.

As of now, default version computes and displays all tidal displacements and stress with Maxwell rheology, and calculates failure using the Mohr Coulomb failure criteria assuming seismic conventions and lithostatic pressure.
If changes are needed, please contact corresponding author at laurent.pou@jpl.nasa.gov

Example files from Pou et al. 2023 - Tidal seismicity in the silicate interior of Europa are given for reference (Europa...txt and Europa_e, and Moon...txt and Moon_e)
Typical run for Europa models with single cohesion and friction values is between 10 to 30 minutes.

Main code functions are stored in the Matlab files folder.
Visualisation tools are stored in both the Matlab files folder and the private subfolder.


Notes:
- Error messages for imaginary parts ignored are only related to graphics and have no impact on visco-elastic calculations.
- Figure 7 will freeze and stay blank until the stress calculations (longest part of the code) are done. 
This will be corrected in a later version.
- Stress calculations are really expensive in terms of memory storage (more than 20 4D matrix are used). 
Be careful and use low sampling values for initial tests.
Memory and Time optimization will be done in future releases.
- There are many figures in the default version (almost 100).
A separate file for figure management will be done soon (currently managed with booleans in Matlab files/Main_Viagra.m).
- This was coded in Matlab2022b, but should work with Matlab 2019a or more recent. Please contact the corresponding author in case of backward compatibility issues.


Corresponding author for the current VIAGRA version and maintenace: Laurent Pou (laurent.pou@jpl.nasa.gov)

Credits for VIAGRA main code: Nicolas Guyennon, Ozgur Karatekin, Laurent Pou

Credits for visualisation tools: M_Map - mapping toolbox (Author: rich@eos.ubc.ca), Version 1.4h  Nov 2014

A paper is in progress for future citation purposes.
