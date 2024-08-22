

K_b = linspace(0.0002,0.00095,20); %range of control parameter defined
tau_b = linspace(0.45,0.85,20); %range of time delay defined

bifurcation(K_b,tau_b)

K_b = [0.0002,0.0007];
X_i = [0.15 0.2];
vel_energy(K_b,0.45,X_i) to python