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
%BC.Boundary.Values = [1 1 1 1];
BC.Boundary.Values = [1,2,4,3];

BC.InputVertexValues = [5 7 0 9];
%BC.InputVertexValues = [3 5 7 9];
% Questi indici posso essere anche indici ai valori numerici
% contenuti nel vettore BC.Values

BC.Holes.Hole = [];
BC.Segments.Segment = [];

mu=1;
beta=[0;0];
gamma=0;
f=@(t,x) 1;

boundaryFunctions{1}= @(t,x) zeros(size(x(1,:)));
boundaryFunctions{2}= @(t,x) zeros(size(x(1,:)));
boundaryFunctions{4}= @(t,x) 10*t*ones(size(x(1,:)));
boundaryFunctions{3}= @(t,x) zeros(size(x(1,:)));

initialSolution=@(x) zeros(size(x(1,:)));
%u=  @(x) 16*x(1,:).*(1-x(1,:)).*x(2,:).*(1-x(2,:));
%grad_u =  @(x)[16*(1-2*x(1,:)).*x(2,:).*(1-x(2,:));16*x(1,:).*(1-x(1,:)).*(1-2*x(2,:))];

problem=ParabolicProblem(mu,beta,gamma,f,...
                         Domain,[0,1],BC,boundaryFunctions,initialSolution);%,...
                                                %u,grad_u);



