[y sllreturn bwfn]= antennafunccircular(x1, null, phi_desired, distance)
is the main function.

Input Description:

x1=  input string 
null= the angles at which the sidelobe level is required to be null controlled.
phi_desired= the desired angle of maximum radiation
distance= Spacing between elements as a multiple of wavelength
Output Description:
y= Objective function returned that needs to be minimized
sllreturn= the maximum sidelobe level is returned
bwfn= the beamwidth between first nulls is also returned

Encoding of input string x1

If the number of elements= N
Dimension of input string= N
The first N/2 dimensions encode the element excitations with bounds (1,N) and the next N/2 dimensions encode the phase excitations with bounds (-180,180).
Note: Here we have not considered the Dynamic Range Ratio as a relevant factor.

display_plot(gbest,phi_desired,distance)  gives us the radiation pattern plot


Sample Inputs:

x1= Any string within bounds
null= [50, 120] in radians
phi_desired= 180 degrees
distance= 0.5


