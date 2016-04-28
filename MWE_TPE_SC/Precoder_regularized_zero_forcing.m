function G_rzf=Precoder_regularized_zero_forcing(P,H,csi)
%The definitions and equations refer (mostly) to the 
%following publication:
%
%Müller, A., A. Kammoun, E. Björnson, and M. Debbah, "Linear Precoding 
%Based on Polynomial Expansion: Reducing Complexity in Massive MIMO", 
%EURASIP Journal on Wireless Communications and Networking, 2016(1), 1-22, 
%DOI: 10.1186/s13638-016-0546-z.
%Pre-print available at: http://www.laneas.com/axel-muller
%
%This is version 1.0. (Last edited: 2016-03-16)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

[M,K]=size(H);
% Q=inv(H*H'+csi*eye(M));
% beta=sqrt(P/(1/K*trace(H'*Q^2*H)));
% G_rzf=beta*inv(H*H'+csi*eye(M))*H;

G_rzf_tmp=inv(H*H'+csi*eye(M))*H;
beta=sqrt(P/(1/K*trace(G_rzf_tmp*G_rzf_tmp')));
G_rzf = beta*G_rzf_tmp;