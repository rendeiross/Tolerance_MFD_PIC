# Tolerance_MFD_PIC


The equations used were retrieved from an open source paper: https://doi.org/10.1364/AO.540682, published at Optica Group 2024.
--------------------////////////////////////-----------------------------------
Running

Run the script and wait for the GUI to pop.

Input the values of Dx1,Dy1 and Dx2,Dy2, representing the diameter of the beam for x and y transversal plane of the coupling. For symmetric beam (circular optical port such as fiber, check symmetric) D1 is port 1, D2 is port 2.

Input the tolerancing values. The wavelength is required for ztolerance calculations.

--------------------////////////////////////-----------------------------------
Results explanation

η_M represents only the insertion loss due to mode mismatch

η_δx and η_δy represent the insertion loss due to transversal misalignment

η_ZM is the η_M + the tolerance in z (optical axis). to obtain η_δz = η_ZM - η_M.

η_total = sums the η_δx, η_δy and η_ZM insertion loss


--------------------////////////////////////-----------------------------------

Generate IL map (2D): generates the heat map from the coupling, drawing the 1 dB, 2dB and 3 dB loss areas.

Generate IL map (1D): same as Generate IL map (2D) but with separated axes to see the difference between x and y.

Generate IL vs deltaz: generates a plot with IL vs ztolerance for a givem max z and step.

Labels: allows to change the names of the labels for easy and fast figure creation. (Generate IL vs deltaz is not integrated in labels, therefore, to change labels for these graphs you need to manually change in the script and re-run the app).


