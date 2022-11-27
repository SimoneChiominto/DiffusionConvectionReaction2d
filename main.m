clear all
%if(~exist('bbtr30'))
%addpath('../triangolatore/Long/bbtr30')
%addpath('../triangolatore/Long')
%addpath("@Mesh")
%addpath("@DiffusionConvectionReactionProblem2D")
disp('mesher added to the path')
%end

AreaValue=0.1;


% -------------------------------
% Inserimento dei vertici
% -------------------------------

Domain.InputVertex = [ 0 0
                       1 0
                       1 1
                       0 1];



% ---------------------------------------------
% Definizione del dominio a partire dai Vertici
% ---------------------------------------------

% Dichiaro le variabili per delimitare il dominio
Domain.Boundary.Values = 1:4;
% lato di bordo 1 dal nodo 1 al nodo 2
% lato di bordo 2 dal nodo 2 al nodo 3
% lato di bordo 3 dal nodo 3 al nodo 4
% lato di bordo 4 dal nodo 4 al nodo 1

Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio

% --------------------------------------------------
% Definizione delle condizioni al contorno a partire
% dai Vertici e dai lati di bordo
% --------------------------------------------------

% valori numerici per le condizioni al contorno
BC.Values = [0.0 2.0 0.0 40.0 0.0 16.0 0.0];

% marker delle condizioni al contorno sui bordi del dominio
% dispari -> Dirichlet; pari -> Neumann

%BC.Boundary.Values = [1 2 4 3];
BC.Boundary.Values = [1 1 1 1];

% marker dei Vertici iniziali

%BC.InputVertexValues = [5 7 0 9];
BC.InputVertexValues = [9 11 13 15];
% Questi indici posso essere anche indici ai valori numerici
% contenuti nel vettore BC.Values

BC.Holes.Hole = [];
BC.Segments.Segment = [];


% --------------------------------------------
% creazione dell'istanza del problema differenziale
% --------------------------------------------
mu=@(x) ones(size(x(1,:)));
beta=@(x) zeros(2,length(x(1,:)));
gamma=@(x) zeros(size(x(1,:)));
f= @(x) 32*(x(1,:).*(1-x(1,:)) + x(2,:).*(1-x(2,:))); % forcing term

boundaryFunctions{1} = @(x) zeros(size(x(1,:)));
boundaryFunctions{2} = @(x) zeros(size(x(1,:)));
boundaryFunctions{4} = @(x) zeros(size(x(1,:)));
boundaryFunctions{3} = @(x) zeros(size(x(1,:)));

u= @(x) 16*x(1,:).*(1-x(1,:)).*x(2,:).*(1-x(2,:)); 
grad_u = @(x)[16*(1-2*x(1,:)).*x(2,:).*(1-x(2,:)); 16*x(1,:).*(1-x(1,:)).*(1-2*x(2,:))];

problem=DiffusionConvectionReactionProblem2D(mu,beta,gamma,f,...
                                                Domain,BC,boundaryFunctions,...
                                                u,grad_u);



% --------------------------------------------
% Inserimento dei parametri di triangolazione
% --------------------------------------------

RefiningOptions.CheckArea  = 'Y';
RefiningOptions.CheckAngle = 'Y';
RefiningOptions.AreaValue  = AreaValue; %l'ho modificato
RefiningOptions.AngleValue = 20;
RefiningOptions.Subregions = [];

% --------------------------------------------
% Creazione della triangolazione
% --------------------------------------------
mesh=Mesh(Domain,BC,RefiningOptions);

% --------------------------------------------
% creazione del elemento di riferimento (Courant)
% --------------------------------------------

courantEl.phi={@(x) x(1,:), @(x) x(2,:),@(x) 1-x(1,:)-x(2,:)};
courantEl.gradPhi={@(x) [1;0], @(x) [0;1],@(x) [-1;-1]};
N1=@(x) x(1,:);
N2=@(x) x(2,:);
N3=@(x) 1-x(1,:)-x(2,:);
grad_N1=[1;0];
grad_N2=[0;1];
grad_N3=[-1;-1];

P2El.phi={@(x) 2*N1(x).*(N1(x)-0.5),...
          @(x) 2*N2(x).*(N2(x)-0.5),...
          @(x) 2*N3(x).*(N3(x)-0.5),...
          @(x) 4*N1(x).*N3(x),...
          @(x) 4*N2(x).*N1(x),...
          @(x) 4*N2(x).*N3(x)...
          };

P2El.gradPhi={@(x) grad_N1*(4*N1(x)-1),...
              @(x) grad_N2*(4*N2(x)-1),...
              @(x) grad_N3*(4*N3(x)-1),...
              @(x) 4*(grad_N1*N3(x)+grad_N3*N1(x)),...
              @(x) 4*(grad_N2*N1(x)+grad_N1*N2(x)),...
              @(x) 4*(grad_N2*N3(x)+grad_N3*N2(x))...
              };


% --------------------------------------------
% creazione dell'istanza del problema approssimato
% --------------------------------------------
approx_problem=ApproxDiffusionConvectionReactionProblem2D(problem,mesh,P2El);


% --------------------------------------------
% risoluzione del problema approssimato
% --------------------------------------------
approx_problem.generateLinearSystem();
approx_problem.solve()
approx_problem.plot()

% --------------------------------------------
% calcolo dell'errore
% --------------------------------------------
approx_problem.getL2Error()
approx_problem.getH0Error()
approx_problem.getLInfError()


% -------------------------------------------
% verifico ordine di convergenza dell'errore
% -------------------------------------------


AreaValue=[0.05 0.03 0.02 0.01 0.008 0.007 0.005 0.003 0.002];% 0.001 0.0005];% 0.002];
max_k=length(AreaValue);

for k=1:max_k

    RefiningOptions.AreaValue  = AreaValue(k);
    mesh=Mesh(problem.domain,problem.BC,RefiningOptions);
    k
    approx_problems(k)=ApproxDiffusionConvectionReactionProblem2D(problem,mesh,courantEl);
    tic
    approx_problems(k).generateLinearSystem();
    toc

    tic
    approx_problems(k).solve()
    toc
    approx_problems(k).plot()
    
    L2_error(k)=approx_problems(k).getL2Error();
    H0_error(k)=approx_problems(k).getH0Error();
    %Linf_error(k)=approx_problems(k).getLInfError();

    h(k)=approx_problems(k).mesh.getDiamMax();
    Area(k)=approx_problems(k).mesh.getAreaMax();
    n_dof(k)=approx_problems(k).mesh.getDoF();

%AreaValue=AreaValue/2;
end

polyfit(log(h),log(L2_error),1)
figure 
loglog(h,L2_error,"o")
%polyfit(log(h),log(Linf_error),1)
%figure 
%loglog(h,Linf_error,"o")
polyfit(log(h),log(H0_error),1)
figure 
loglog(h,H0_error,"o")

polyfit(log(Area),log(L2_error),1)
figure 
loglog(Area,L2_error)
polyfit(log(Area),log(H0_error),1)
figure 
loglog(Area,H0_error)

polyfit(log(n_dof),log(L2_error),1)
figure 
loglog(n_dof,L2_error)
polyfit(log(n_dof),log(H0_error),1)
figure 
loglog(n_dof,H0_error)




