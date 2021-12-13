# MultiMachineDAEControl
The codes correspond to our paper on load frequency control for power systems DAE which is accessible here https://arxiv.org/abs/2104.05957. The main program is the main code to run the simulation on the 3-machine, 9-bus system. To run it, you'll need MATPOWER to solve the power flow problem, PST at least for the description of generator parameters, and YALMIP with MOSEK solver to solve LMIs.
