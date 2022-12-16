clear
Domain.InputVertex = [ 0 0
                       1 0
                       1 1
                       0 1];
Domain.Boundary.Values = 1:4;
Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio

BC.Values = [0.0 2.0 0.0 40.0 0.0 16.0 0.0];
% marker delle condizioni al contorno sui bordi del dominio
% dispari -> Dirichlet; pari -> Neumann
BC.Boundary.Values = [1 1 1 1];

%BC.InputVertexValues = [5 7 0 9];
BC.InputVertexValues = [3 5 7 9];
% Questi indici posso essere anche indici ai valori numerici
% contenuti nel vettore BC.Values

BC.Holes.Hole = [];
BC.Segments.Segment = [];

mu=1;%@(x) 0.0001%*ones(size(x(1,:)));
%mu=0.000001;
%mu=0.00001;
%beta=0*[100000;0];%@(x)[cos(-0.5+x(1));sin(-0.5+x(2))]/norm(-0.5+x)%*ones(2,length(x(1,:)));
%beta=@(x) 1000000*[sin(x(1)*x(2));x(1)*exp(-x(2))];
beta=[0;0];
%beta=@(x) 1000000*[sin(x(1,:).*x(2,:));exp(x(2,:))];%1000000*[1;1];
gamma=0;
%gamma=0*1000000;%@(x) 1* ones(size(x(1,:)));
%gamma=@(x) 1000000*log(3+(x(1)*x(2))^2);
%f=@(x) 0.0001*4*pi^2*(x(1,:).^2+x(2,:).^2).*sin(2*pi*x(1,:).*x(2,:))...
%    +1000000*2*pi*(sin(x(1,:).*x(2,:)).*x(2,:)+exp(x(2,:)).*x(1,:))*cos(2*pi*x(1,:).*x(2,:));
    %+ 1000000*sin(2*pi*x(1,:).*x(2,:));
%f= @(x) 0.0001*(x(1,:).^2+x(2,:).^2).*sin(x(1,:).*x(2,:))+...
%    1000000*(x(1,:)+x(2,:)).*cos(x(1,:).*x(2,:));
f= @(x) 32*(x(1,:).*(1-x(1,:)) + x(2,:).*(1-x(2,:))); % forcing term
%f=@(x) -0.000018*x(1)*x(2)+90000*x(1)^2*x(2)-(10000+30000*x(1)^3)
%f=@(x) -0.00018*x(1)*x(2)+1000000*(9*x(1)^2*x(2)*sin(x(1)*x(2))+x(1)*exp(-(x(2)))*(1+3*x(1)^3))
%f=@(x) -0.00018*x(1)*x(2)+1000000*log(3+(x(1)*x(2))^2)*(1+3*x(1,:).^3)*x(2,:)



%boundaryFunctions{1} = @(x) sin(2*pi*x(1,:).*x(2,:));%
boundaryFunctions{1}= @(x) x(1,:).*x(2,:);
%boundaryFunctions{1} = @(x) zeros(size(x(1,:)));
%boundaryFunctions{2} = @(x) zeros(size(x(1,:)));
%boundaryFunctions{4} = @(x) zeros(size(x(1,:)));
%boundaryFunctions{3} = @(x) x(1,:)+ x(2,:).^2+zeros(size(x(1,:)));

%u= @(x) sin(2*pi*x(1,:).*x(2,:));%
u= @(x) 16*x(1,:).*(1-x(1,:)).*x(2,:).*(1-x(2,:))+x(1,:).*x(2,:); 
grad_u = @(x)[16*(1-2*x(1,:)).*x(2,:).*(1-x(2,:))+x(2,:);...
              16*x(1,:).*(1-x(1,:)).*(1-2*x(2,:))+x(1,:)];
%grad_u = @(x)[2*pi*x(2,:).*cos(2*pi*x(1,:).*x(2,:))+x(2,:);...
%              2*pi*x(1,:).*cos(2*pi*x(1,:).*x(2,:))+x(1,:)];
problem2=DiffusionConvectionReactionProblem2D(mu,beta,gamma,f,...
                                                Domain,BC,boundaryFunctions,...
                                                u,grad_u);

