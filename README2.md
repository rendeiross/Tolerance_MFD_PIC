# Tolerance_MFD_PIC
This README is for Optical coupling calculator_v2.py

The equations used were retrieved from an open source paper: https://doi.org/10.1364/AO.540682, published at Optica Group 2024.


## Running

Run the script and wait for the GUI to pop.


## Tab 1: Coupling Calculator

Input the values of Dx1,Dy1 and Dx2,Dy2, representing the diameter of the beam for x and y transversal plane of the coupling. For symmetric beam (circular optical port such as fiber, check symmetric) D1 is port 1, D2 is port 2.

Input the tolerancing values. The wavelength is required for ztolerance calculations.

Input n1 and n2 correspond to refractive index of optical port 1 and 2 respectively

Changing the values in input changes the visualization on the right. 

Results explanation

η_M represents only the insertion loss due to mode mismatch

η_δx and η_δy represent the insertion loss due to transversal misalignment

η_ZM is the η_M + the tolerance in z (optical axis). to obtain η_δz = η_ZM - η_M.

η_Fres represents the loss due to fresnel reflections

Total = sums the η_δx, η_δy and η_ZM and fresnel insertion loss

tick dB for results in dB, keep unticked for results in percentage.

tick +Fresnel to add Fresnel losses due to n1 and n2 interface. If z=0 fresnel loss is due to n1/n2 interface. if z>0 fresnel loss is n1/air + air/n2

Labels: allows to change the names of the labels for easy and fast figure creation.


PLOTS:
Top left plot shows the cross-section of both optical ports and the corresponding misalignment.
Top right plot shows the heatmap of the insertion loss due to mode mismatch, draws lines for -1dB and -3dB tolerance and sets a cross at the chosen deltax,deltay point.
Bottom left plot shows the 1D misalignment vs insertion loss for x and y sweep. Places a point at the chosen deltax,deltay.
Bottom left plot shows the delta z vs insertion loss. Places a point at the chosen deltaz.

## Tab 2: Converter

Simple converter percentage to dB and dB to percentage. current to optical power. dBm/mW conversion and half angle divergence calculations

## Tab 3: Fresnel stack

A simple APP using the formula R = ((n1 - n2) / (n1 + n2)) ** 2, where the user can stack some materials and evaluate the Fresnel losses due to reflections at the irterface. The user can quickly choose the layer in the right menu and use keyboard arrows up and down to move the different materials.
The user can choose the predefined materials or simply write a new one, since the name and refractive index boxes can be interactively changed.

The fresnel_stack.py can be runned 'solo' or in Optical coupling calculator.py, integrated as tab 3 (both need to be on the same folder)

