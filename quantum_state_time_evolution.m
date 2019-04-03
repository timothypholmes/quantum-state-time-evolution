%%
%                                                                                                                                
%                   Timothy Holmes                            
%             quantum_state_time_evolution                                                    
%                      (3/2/18)                               
%                                                                                                         
%% Function

S = load('C-2-Data.mat');
H = A;
psi_0 = state1;

function quantum_state_time_evolution(time,H,psi_0)
%% load sample
S = load('C-2-Data.mat');
%[ 0 -1i 0 ; 1i 0 -1i; 0 1i 0 ]/sqrt(2),[ 3 ; 1i ; -1i ]
%[ 3 ; 1i ; -1i ]
%% Eigen terms

hbar = 1; %hbar set to a value of 1
[V,D] = eig(H); %Used to find the eigenvectors and the eigen 
%values on the diagonal of the identity matrix. This can also be 
%found by the equation D = P*H*P^-1.
eigVal = eig(H); %This finds all the eigen values of H and puts them 
%into a column vector ( Nx1 ).
eStates = V; %renameing the eigenvectors
diag = D; %renaming the diagonal matrix

%% Normalizing incoming state

norms = norm(psi_0); %using the Norm function to the normalized vector.
%This finds magnitide of the vector by taking the inner product. 
normState = (1/norms); %To normalize the state we have to find 1/norm
Normin = normState.*psi_0; %Then take 1/norm and multiply it by the original
%incoming state. 

%% Finding C_n terms

c_n = zeros(1,length(eStates)); %Sets a row vector (1xN) filled with zeros

for j = 1:length(Normin) 
    c_n(j) = dot(eStates(:,j),Normin);
    
end 

%% Finding time evolution psi_t

psi_t = zeros(length(eStates),length(time)); %presets a NxM matrix to be
%filled with zeros.

count = 0; %sets the values of count equal to zero.

for k = 1:length(time) %nested for loop with 
    for i = 1:length(Normin)
 
    psi_t(:,k) = psi_t(:,k) + (c_n(1,i).*(exp((-1i*eigVal(i,1).*time(k))/hbar))).*eStates(:,i);

    count = count + 1; % Counts how many iterations are executed in the for loop

    end 
    
end

%% Plot

figure(1) 
hold on;
    plot(time,real(psi_t(1,:)))
    plot(time,imag(psi_t(1,:)))
legend('real','imaginary')
title('Time Evolution')
xlabel('Time')
ylabel('\bf{\psi(t)}')
fprintf('The number of iterations: %f \n', count)


end





























