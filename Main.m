clear
clc

%I am here introducing the variables

%The number of meshes to be inputted by the user "E" will be used
E= input("hi doctor please input the number of elements you desire ");

tic

%The wall thickness it is 0.4 by default from the problem
L=0.4;

%the number of nodes in the system
N=E+1;

%The number of nodes within every element it will be 4 since we are using
%cubic elements
n=4;

%The length of each element is
he=L/E;

%The number of Gauss points needed for my integral and since my integral
%has an order 5, and using nGP=(order+1)/2
nGP=3;

%I am introducing the K function which is the general coefficient and no
%element wise

%I am here introducing the shape functions
psi_coeff=[-9/(2*(he^3)), 9/(he^2), -11/(2*he), 1;
           27/(2*(he^3)), -45/(2*(he^2)), 9/he, 0;
           -27/(2*(he^3)), 18/(he^2), -9/(2*he), 0;
           9/(2*(he^3)), -9/(2*(he^2)), 1/he, 0];

%I am here introducing the derivatives of each chape function
psi_dot_coeff=[-27/(2*(he^3)), 18/(he^2), -11/(2*he);
                81/(2*(he^3)), -45/(he^2), 9/he;
               -81/(2*(he^3)), 36/(he^2), -9/(2*he);
                27/(2*(he^3)), -9/(he^2), 1/he];
           
%I am defining the Gauss points
Gauss_points=[0, 0.8888888889;
             -0.7745966692, 0.5555555555;
              0.7745966692, 0.5555555555];
          
%Creating the B matrix or connectivity matrix
connectivity=1;
B=zeros(E,n);
for r=1:E
    for s=1:n-1:n
        B(r,s)=connectivity;
        connectivity=connectivity+1;
    end
    connectivity=connectivity-1;
end
connectivity=connectivity+1;
for r=1:E
    for s=2:1:n-1
        B(r,s)=connectivity;
        connectivity=connectivity+1;
    end
end

%The number of total nodes within the entire system including virtual nodes
NN=max(max(B));

%I am introducing the K function which is the general coefficient and not
%element wise
K=sparse(NN,NN);

%I am now initializing xa
xa=0;

for e=1:E
    for i=1:n
        for j=1:n
            k_for_iteration=0;
            for GP=1:nGP
                
                %Mapping between between x_bar and zeta, and zeta is all
                %the points at the first column of the Gauss_points matrix
                x_bar=(Gauss_points(GP,1)+1)*he/2;
                
                %Computing the thermal conductivity
                thermal_conductivity=10*(x_bar+xa)+1;
                
                %Computing the actual psi to work with
                psi_dot_coeffs_p2=[x_bar^2; x_bar; 1];
                psi_dot_j=(psi_dot_coeff(j, :)*psi_dot_coeffs_p2);
                psi_dot_i=(psi_dot_coeff(i, :)*psi_dot_coeffs_p2);
                
                %Calculating the Gauss point weight the weight is the
                %second column of the Gauss_points matrix
                
                %Computing k for a single Gauss point and adding it to the
                %previous one to compute the integral
                k=-1*thermal_conductivity*psi_dot_j*psi_dot_i*(he/2)*Gauss_points(GP,2);
                k=k_for_iteration+k;
                k_for_iteration=k;
            end
            K(B(e,i),B(e,j))=k+K(B(e,i),B(e,j));
            
        end
    end
    xa=xa+he;
end

%I am now working on the boundary conditions I will create the Q matrix in
%due time, I want to be working on one thing at a time

%step one to introduce the first dericle boundary condition I will zero the
%entire first line, and make sure that the K11 is equal 1
K(1,:)=0;
K(1,1)=1;

%I will now be working on the Q matrix, the left side really, and since F
%is equal to 0 I am happy, no worries all is easy it is simply all equal to
%0, with me only having to play with the boundary conditions

Q=zeros(NN, 1);
Q(1)= 15;
Q(N)=10;

%Now I am introducing the solutions under the form U
U_pre_processing=K\Q;

%I am introducing the U_per_element_prime matrix THE U_per_element_prime IS
%NOT A DERIVATIVE IT IS ONLY A MATRIX TO SIMPLIFY CALCULATIONS
U_per_element=zeros(E,n);

%I am creating the U_per_element_prime matrix
for e=1:E
    for i=1:n
U_per_element(e,i)=U_pre_processing(B(e,i));
    end
end

%I am finding the number of nodes in the port processing
post_processing_nodes_per_element=40;
post_processing_nodes=post_processing_nodes_per_element*E-(E-1);

%I am creating the x bar matrix
x_bar=[(linspace(0, he, post_processing_nodes_per_element).^3);
       (linspace(0, he, post_processing_nodes_per_element).^2);
       (linspace(0, he, post_processing_nodes_per_element));
       (linspace(1, 1, post_processing_nodes_per_element))];

%I am creating the psi functions         
psi_functions=psi_coeff*x_bar;

%Initialization
U=zeros(post_processing_nodes, 1);

%Now I am creating the U matrix
U_prime=U_per_element*psi_functions;
U(1:post_processing_nodes_per_element)=U_prime(1,:);
element_start=post_processing_nodes_per_element+1;
element_end=post_processing_nodes_per_element*2-1;
for element=2:E
    U(element_start:1:element_end)=U_prime(element,2:1:post_processing_nodes_per_element);
    element_start=element_start+post_processing_nodes_per_element-1;
    element_end=element_end+post_processing_nodes_per_element-1;
end

%entering U and plotting
U=[transpose(linspace(0,L,post_processing_nodes)),U];
plot(U(:,1),U(:,2))

%plot design
xlabel_me=xlabel('x (m)');
ylabel_me=ylabel('Temperature (C)');

A1=toc;

%new part

Uold=U;

%I am here introducing the variables

%The wall thickness it is 0.4 by default from the problem
L=0.4;

%the number of nodes in the system
N=E+1;

%The number of nodes within every element it will be 4 since we are using
%cubic elements
n=4;

%The length of each element is
he=L/E;

%The number of Gauss points needed for my integral and since my integral
%has an order 5, and using nGP=(order+1)/2
nGP=3;

%I am introducing the K function which is the general coefficient and no
%element wise

%I am here introducing the shape functions
psi_coeff=[-9/(2*(he^3)), 9/(he^2), -11/(2*he), 1;
           27/(2*(he^3)), -45/(2*(he^2)), 9/he, 0;
           -27/(2*(he^3)), 18/(he^2), -9/(2*he), 0;
           9/(2*(he^3)), -9/(2*(he^2)), 1/he, 0];

%I am here introducing the derivatives of each chape function
psi_dot_coeff=[-27/(2*(he^3)), 18/(he^2), -11/(2*he);
                81/(2*(he^3)), -45/(he^2), 9/he;
               -81/(2*(he^3)), 36/(he^2), -9/(2*he);
                27/(2*(he^3)), -9/(he^2), 1/he];
           
%I am defining the Gauss points
Gauss_points=[0, 0.8888888889;
             -0.7745966692, 0.5555555555;
              0.7745966692, 0.5555555555];
          
%Creating the B matrix or connectivity matrix
connectivity=1;
B=zeros(E,n);
for r=1:E
    for s=1:n-1:n
        B(r,s)=connectivity;
        connectivity=connectivity+1;
    end
    connectivity=connectivity-1;
end
connectivity=connectivity+1;
for r=1:E
    for s=2:1:n-1
        B(r,s)=connectivity;
        connectivity=connectivity+1;
    end
end

Nzeroes=(n*n*E-NN);

%The number of total nodes within the entire system including virtual nodes
NN=max(max(B));

%I am introducing the K function which is the general coefficient and not
%element wise
Knew=sparse(NN,NN, Nzeroes);

%I am now initializing xa
xa=0;

%Mapping between between x_bar and zeta, and zeta is all the points at the
%first column of the Gauss_points matrix
                x_bar_a=(Gauss_points(:,1)+1).*he/2;

%Computing the actual psi to work with
                psi_dot_coeffs_p2=[(x_bar_a').^2; (x_bar_a'); (x_bar_a')./(x_bar_a')];
                psi_dot=psi_dot_coeff*psi_dot_coeffs_p2;

for e=1:E
    for i=1:n
        for j=i:n
            k_for_iteration=0;
            for GP=1:nGP
                x_bar=(Gauss_points(GP,1)+1).*he/2;

                
                %Computing the thermal conductivity
                thermal_conductivity=10*(x_bar+xa)+1;
                
             
                
                %Calculating the Gauss point weight the weight is the
                %second column of the Gauss_points matrix
                
                %Computing k for a single Gauss point and adding it to the
                %previous one to compute the integral
                Knew(B(e,i),B(e,j))=-1*thermal_conductivity*psi_dot(i,GP)*psi_dot(j,GP)*(he/2)*Gauss_points(GP,2)+Knew(B(e,i),B(e,j));
            end
            
        end
    end
    xa=xa+he;
end



Knew = (Knew+Knew') - eye(size(Knew)).*Knew;

%I am now working on the boundary conditions I will create the Q matrix in
%due time, I want to be working on one thing at a time

%step one to introduce the first dericle boundary condition I will zero the
%entire first line, and make sure that the K11 is equal 1
Knew(1,:)=0;
Knew(1,1)=1;

%I will now be working on the Q matrix, the left side really, and since F
%is equal to 0 I am happy, no worries all is easy it is simply all equal to
%0, with me only having to play with the boundary conditions

Qnew=zeros(NN, 1);
Qnew(1)= 15;
Qnew(N)=10;

%Now I am introducing the solutions under the form U
U_pre_processing=Knew\Qnew;

%I am introducing the U_per_element_prime matrix THE U_per_element_prime IS
%NOT A DERIVATIVE IT IS ONLY A MATRIX TO SIMPLIFY CALCULATIONS
U_per_element=zeros(E,n);

%I am creating the U_per_element_prime matrix
for e=1:E
    for i=1:n
U_per_element(e,i)=U_pre_processing(B(e,i));
    end
end

%I am finding the number of nodes in the port processing
post_processing_nodes_per_element=40;
post_processing_nodes=post_processing_nodes_per_element*E-(E-1);

%I am creating the x bar matrix
x_bar=[(linspace(0, he, post_processing_nodes_per_element).^3);
       (linspace(0, he, post_processing_nodes_per_element).^2);
       (linspace(0, he, post_processing_nodes_per_element));
       (linspace(1, 1, post_processing_nodes_per_element))];

%I am creating the psi functions         
psi_functions=psi_coeff*x_bar;

%Initialization
U=zeros(post_processing_nodes, 1);

%Now I am creating the U matrix
U_prime=U_per_element*psi_functions;
U(1:post_processing_nodes_per_element)=U_prime(1,:);
element_start=post_processing_nodes_per_element+1;
element_end=post_processing_nodes_per_element*2-1;
for element=2:E
    U(element_start:1:element_end)=U_prime(element,2:1:post_processing_nodes_per_element);
    element_start=element_start+post_processing_nodes_per_element-1;
    element_end=element_end+post_processing_nodes_per_element-1;
end

%entering U and plotting
U=[transpose(linspace(0,L,post_processing_nodes)),U];
plot(U(:,1),U(:,2))

%plot design
xlabel_me=xlabel('x (m)');
ylabel_me=ylabel('Temperature (C)');

A2=toc;

A3=(A1-(A2-A1))/A1;
disp(A3);

disp(sqrt(sum((U(:,2)-Uold(:,2)).^2)/sum(U(:,2).^2)));