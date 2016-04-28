function y=fixed_point_solve(Phi,tau,rho_lin,K)
%The references to definitions and equations refer (mostly) to the 
%following publication:
%
%Müller, A., A. Kammoun, E. Björnson, and M. Debbah, "Efficient Linear 
%Precoding for Massive MIMO Systems using Truncated Polynomial Expansion", 
%IEEE 8th Sensor Array and Multichannel Signal Processing Workshop (SAM), 
%Coruna, Spain, 2014.
%Available at: http://www.laneas.com/axel-muller
%
%This is version 1.1. (Last edited: 2016-03-14)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

csi=0.5;
y=equation_solve(Phi,csi,tau,rho_lin,K);


while(abs(y-csi)/csi>1e-4)  % Check for convergence; 1e-4 should be enough
csi=y;
y=equation_solve(Phi,csi,tau,rho_lin,K);
end

