# NonLocalTunnelingModel
This is the Non-Local-Tunneling Model used for simulating source to drain tunneling in ultra short double-gate transistors.

1. M0.mat-M11.mat or M0.mat-M11.csv for future python implementations.
The mat files are postprocessing files from Sentaurus S-Device from a double gate InGaAs ultra thin film FET, Lg=11.5 nm.
The pre-processing files include the conduction band (CB) energies, the fermi energies(FE), Quasi-Fermi Potentials (QSFP) and the current (Ic).

2. tunneling2finnerQSFP.m 
Calculates the Generatio and Recombination Rate for the Non-Local Tunneling Model, depending the QSFP (line 108-109) at teach point of the tunneling path or at the Source/Drain Contact.



