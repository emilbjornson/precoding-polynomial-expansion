function [gf,sinr,A_bar,B_bar,C_bar]=Precoder_polynomial_updated(P,Phi,L,tau,sigma,K)
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

M=length(Phi);
derivation_order=2*(L-1)+1;
delta_tab=zeros(derivation_order+1,1);
delta_tab(1)=1/K*trace(Phi);
f_tab=zeros(derivation_order+1,1);
g_tab=zeros(derivation_order+1,1);
j_tab=zeros(derivation_order+1,1);
g_tab(1)=0;
f_tab(1)=-1/(1+g_tab(1));
T=zeros(M,M,derivation_order+1);
Q=zeros(M,M,derivation_order+1);
T(:,:,1)=eye(M);
Q(:,:,1)=zeros(M,M);
%f_tab(1)*Phi;
for k=1:derivation_order

	Q(:,:,k+1)=k*f_tab(k)*Phi;
	for i=0:k-1
		for j=0:i
		T(:,:,k+1) = T(:,:,k+1) + nchoosek(k-1,i)*nchoosek(i,j)*T(:,:,k-i)*Q(:,:,i-j+2)*T(:,:,j+1);
		end
	end
	scalar=0;
	for i=0:k-1
		for j=0:i
			scalar=scalar+nchoosek(k-1,i)*nchoosek(i,j)*(k-i)*f_tab(j+1)*f_tab(i-j+1)*delta_tab(k-i);
		end
	end
	f_tab(k+1)=scalar;
	g_tab(k+1)=k*delta_tab(k);
	delta_tab(k+1)=1/K*trace(Phi*T(:,:,k+1));
end

j_tab(1)=delta_tab(1)*f_tab(1);
for k=1:derivation_order
	for n=0:k
		j_tab(k+1)=j_tab(k+1)+nchoosek(k,n)*delta_tab(n+1)*f_tab(k-n+1);
	end
end
gj=zeros(derivation_order+1,1);
gj(1)=g_tab(1)*j_tab(1);
for k=1:derivation_order
	for n=0:k
		gj(k+1)=gj(k+1)+nchoosek(k,n)*g_tab(n+1)*j_tab(k-n+1);
	end
end
Z_tab=-sqrt(1-tau^2)*j_tab;
X_tab(1)=delta_tab(1)*1/K*trace(Phi)+1/K*trace(Phi*Phi)-tau^2*delta_tab(1)^2;   % In the paper we have Phi*Phi*T, but T=eye() @t=0
X_tab(2:derivation_order)=1./[2:1:derivation_order]'.*(j_tab(3:end)+tau^2*gj(3:end));

A_tab=(((-1).^[0:L-1]')./factorial([0:L-1]')).*Z_tab(1:L);
A_bar=nan(L);
for l=0:L-1
	for m=0:L-1
		A_bar(l+1,m+1)=A_tab(l+1)*A_tab(m+1);
	end
end
B_bar=nan(L);
for l=0:L-1
	for m=0:L-1
		B_bar(l+1,m+1)=(-1)^(l+m)*X_tab(l+m+1)/factorial(l+m);
	end
end
C_bar=nan(L);
for l=0:L-1
	for m=0:L-1
		C_bar(l+1,m+1)=(-1)^(l+m+1)*1/K*trace(T(:,:,l+m+2))/(factorial(l+m+1));
	end
end


%Updated implementation following Emil
Csqrt = sqrtm(C_bar);

Numerator = (Csqrt\A_bar)/Csqrt;
Denominator =  (Csqrt\(B_bar-A_bar))/Csqrt +(sigma/P)*eye(L);

[g_sub,~] = eigs(Numerator,Denominator,1);
g_sub = g_sub/norm(g_sub); % Circumvent Matlab bug... eigs eigenvectors are not always normalized; causes problems for L=2

gf = sqrt(P)*(Csqrt\g_sub);

sinr=gf'*A_bar*gf/(gf'*(B_bar-A_bar)*gf+sigma);


