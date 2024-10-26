# Thesis_Code

Code attached to the work for the Master Thesis "A convolution process for spatial statistical models on directed linear networks".

All plots and simulations come from these codes.

## Code
### General functions
The two files present relates to the general functions necessary to run all the simulations and the constructions developed.
These are ad hoc created functions, given that there is the need to manipulate the data and the covariance model in an innovative way.
_Functions_Network_ relates to the functions aiming at constructing the network. From the dataframe df, these functions create the actual linear network. Moreover, they create the **paths** object, where all the relevant paths together with the necessary weights are stored.
_Functions_Simulation_ aims to make the actual simulations. In particular, simulate a process, creating a covariance from a paths object. Moreover, the estimation is tackled, functions for all types of estimators are defined here. Lastly, functions for prediction at unsampled locations are defined.