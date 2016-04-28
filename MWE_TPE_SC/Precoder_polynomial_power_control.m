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

function [gf,sinr,A_bar,B_bar,C_bar]=Precoder_polynomial_power_control(Power,P,Phi,L,tau,sigma,K)
% Returns deterministic equivalents for matrices and weights used in power
% controlled TPE.
% Inputs: 
% Power...  Total power of the precoding Matrix.
% P...      Diagonal matrix containing the relative powers of users.
% L...      Derivation order / TPE order.
% tau...    Channel quality parameter.
% sigma... 	Noise variance.
% K...      Number of users.
% Ouputs:
% gf...     Weight vector needed to construct TPE precoding matrix.
% sinr...   Deterministic Equivalent of the SINR achieved using TPE with the gf weight vector.
% A_bar...  Deterministic equivalent of matrix A
% B_bar...  Deterministic equivalent of matrix B
% C_bar...  Deterministic equivalent of matrix C

%% Determine system parameters
M=length(Phi);
derivation_order=L;

%% Algorithm 2: Iterative algorithm for computing \matbf T ^(q)
% Follows J. Hoydis, M. Kobayashi, and M. Debbah "Asymptotic Moments for Interference Mitigation in Correlated Fading Channels", IEEE International Symposium on Information Theory (ISIT'11), Saint-Petersburg, Russia, Jul.-Aug. 2011.
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

%% Algorithm 1: Iterative algorithm for the computation of \beta^(l,m)_M or \nu^(l,m)_M
Tau=zeros(M,M,derivation_order+1);
Tau(:,:,1)=T(:,:,1);
for k=1:derivation_order
	for jj=0:k
		Tau(:,:,k+1)=Tau(:,:,k+1)-nchoosek(k,jj)*T(:,:,jj+1)*f_tab(k-jj+1);
	end
end
beta_vec=zeros(derivation_order+1,derivation_order+1);

for k=1:derivation_order+1
beta_vec(k,1)=1/K*trace(Phi*Tau(:,:,k)*Phi*Tau(:,:,1));
end

for m=2:derivation_order+1
	for k=1:derivation_order+1
		if (k==1), beta_vec(k,m)=1/K*trace(Phi*Tau(:,:,k)*Phi*Tau(:,:,m)); 
		else
			for k_index=0:k-2
				for m_index=0:m-2
					beta_vec(k,m)=beta_vec(k,m)+(k_index+1)*nchoosek(k-1,k_index+1)*(m_index+1)*nchoosek(m-1,m_index+1)*beta_vec(k_index+1,m_index+1)*1/K*trace(Phi*Tau(:,:,m-m_index-1)*Phi*Tau(:,:,k-k_index-1));
				end
			end
		beta_vec(k,m)=beta_vec(k,m)+1/K*trace(Phi*Tau(:,:,k)*Phi*Tau(:,:,m)); 
		end
	end
end

%% Calculate bbar^(l,m)_M
b_vec=zeros(derivation_order+1,derivation_order+1);
for k=0:derivation_order
    for m=0:derivation_order
        for index_k=0:k
            for index_m=0:m
                b_vec(k+1,m+1)=b_vec(k+1,m+1)+(1-tau^2)*nchoosek(k,index_k)*nchoosek(m,index_m)*beta_vec(k-index_k+1,m-index_m+1)*f_tab(index_k+1)*f_tab(index_m+1);
            end
        end
    end
end
for k=1:derivation_order+1
    for m=1:derivation_order+1
        b_vec(k,m)=b_vec(k,m)+tau^2*beta_vec(k,m);
    end
end

%% Calculate Xbar^(l,m)_M
X_k_tab=zeros(derivation_order+1,derivation_order+1);
for k=0:derivation_order
    for m=0:derivation_order
        for index_k=0:k
            for index_m=0:m
                X_k_tab(k+1,m+1)=(1-tau^2)*nchoosek(k,index_k)*nchoosek(m,index_m)*delta_tab(index_k+1)*delta_tab(index_m+1)*f_tab(k-index_k+1)*f_tab(m-index_m+1)+X_k_tab(k+1,m+1);
            end
        end
    end
end

%% Calculate c^(l,m)
c_vec=zeros(derivation_order+1,derivation_order+1);
  for k=1:derivation_order
     for m=1:derivation_order
         for index_k=1:k
             for index_m=1:m
                 c_vec(k+1,m+1)=c_vec(k+1,m+1)+index_k*index_m*nchoosek(k,index_k)*nchoosek(m,index_m)*beta_vec(index_k,index_m)*1/K*trace(Phi*Tau(:,:,k-index_k+1)*Tau(:,:,m-index_m+1));
             end
         end
     end
 end
 for k=0:derivation_order
     for m=0:derivation_order
         c_vec(k+1,m+1)=c_vec(k+1,m+1)+1/K*trace(Phi*Tau(:,:,m+1)*Tau(:,:,k+1));
     end
 end

%% Fill A_bar, B_bar, C_bar
A_bar=zeros(derivation_order+1,derivation_order+1);
B_bar=zeros(derivation_order+1,derivation_order+1);
C_bar=zeros(derivation_order+1,derivation_order+1);

for index_l=0:derivation_order
    for index_c=0:derivation_order
        A_bar(index_l+1,index_c+1)=X_k_tab(index_l+1,index_c+1)*(-1)^(index_l+index_c)/(factorial(index_l)*factorial(index_c));
    end
end

for index_l=0:derivation_order
    for index_c=0:derivation_order
            B_bar(index_l+1,index_c+1)=b_vec(index_l+1,index_c+1)*(-1)^(index_l+index_c)/(factorial(index_l)*factorial(index_c));
    end
end

for index_l=0:derivation_order
    for index_c=0:derivation_order
        C_bar(index_l+1,index_c+1)=c_vec(index_l+1,index_c+1)*(-1)^(index_l+index_c)/(factorial(index_l)*factorial(index_c));
    end
end


%% Find optimal weights
D=B_bar+(sigma/Power)*C_bar;
invsqrtmD=inv(sqrtm(D));
[g_sub,d]=eigs(invsqrtmD*A_bar*invsqrtmD,1);
g_sub = g_sub/norm(g_sub);                          % Fix Matlab eigs normalization bug for 2x2 matrices;
alpha=norm(sqrtm(C_bar)*invsqrtmD*g_sub)^2;
gf=sqrt(Power/(alpha*trace(P)))*invsqrtmD*g_sub;  

%% DE of the SINR
% Using the weights "gf" one should achieve the following SINR in the large system limit.
sinr=zeros(K,1);
for k=1:K
        sinr(k)=K*P(k,k)*gf'*A_bar*gf/(trace(P)*gf'*(B_bar)*gf+sigma);
end
% Remark: 
% There is a shortcoming in the publication, where the DE is stated to be based on 
% Equation (35): sinr(k)=K*P(k,k)*lambdamax/(trace(P));
% While it would have been nicer to state that it is based on Equation "(26.5)":
% sinr(k)=K*P(k,k)*gf'*A_bar*gf/(trace(P)*gf'*(B_bar)*gf+sigma);
% HOWEVER, both of these expressions are ultimatley equivalent! (See below).


%% Debug
% lambdamax = d;
% for k=1:K
%         sinr(k)=K*P(k,k)*lambdamax/(trace(P));
% end

