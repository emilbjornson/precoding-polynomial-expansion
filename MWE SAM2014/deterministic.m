function delta=deterministic(Phi,t,K)
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

Dphi=eig(Phi);
delta=1/K*sum(Dphi);
M=length(Dphi);
delta_prev=0;
while(abs(delta-delta_prev)>1e-8)
delta_prev=delta;
delta=1/K*sum(Dphi./(1+t*Dphi/(1+t*delta_prev)));
end


