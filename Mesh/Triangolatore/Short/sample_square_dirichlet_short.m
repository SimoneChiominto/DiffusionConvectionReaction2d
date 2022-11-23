if(~exist('bbtr30'))
     addpath('../bbtr30')
     disp('../bbtr30 added to the path')
end

%----------------------------------------------------------------------------
%
% Triangolazione di un dominio quadrata
% con condizioni di Dirichlet sul bordo
%
%----------------------------------------------------------------------------
%
%  Autore: Stefano Berrone
%  Politecnico di Torino
%
%----------------------------------------------------------------------------

%clc
%clear all

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
BC.Values = [0.0 12.0 0.0 14.0 0.0 16.0 0.0 0.0 0.0];

% marker delle condizioni al contorno sui bordi del dominio
% dispari -> Dirichlet; pari -> Neumann
BC.Boundary.Values = [3 5 7 9];
% marker dei Vertici iniziali
BC.InputVertexValues = [1 1 1 1];
% Questi indici posso essere anche indici ai valori numerici
% contenuti nel vettore BC.Values

BC.Holes.Hole = [];
BC.Segments.Segment = [];



% --------------------------------------------
% Inserimento dei parametri di triangolazione
% --------------------------------------------

RefiningOptions.CheckArea  = 'Y';
RefiningOptions.CheckAngle = 'N';
RefiningOptions.AreaValue  = 0.02;
RefiningOptions.AngleValue = [];
RefiningOptions.Subregions = [];


% --------------------------------------------
% Creazione della triangolazione e plottaggio
% --------------------------------------------

[geom] = mybbtr30(Domain,BC,RefiningOptions);
draw_grid (geom,1);

% --------------------------------------------------
% --------------------------------------------------


% --------------------------------------------------
% Rielaborazione dei prodotti del triangolatore
% per un piu` agevole trattamento delle condizioni
% al contorno
% --------------------------------------------------

geom.obj.P = geom.obj.P(1:geom.Nobj.N_node,:);
geom.obj.T = geom.obj.T(1:geom.Nobj.N_ele,:);
geom.obj.E = geom.obj.E(1:geom.Nobj.N_edge,:);
geom.obj.Neigh = geom.obj.Neigh(1:geom.Nobj.N_ele,:);

% --------------------------------------------------

j  = 1;
Dj = 1;
for i=1:size(geom.piv.nlist)
     if geom.piv.nlist(i)==0
        geom.piv.piv(i)=j;
        j = j+1;
     else
        geom.piv.piv(i)=-Dj;
        Dj = Dj + 1;
     end
end

% --------------------------------------------------

geom.piv.piv = transpose(geom.piv.piv);

% --------------------------------------------------

% geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
% di Dirichlet e il corrispondente marker

[X,I] = sort(geom.piv.Di(:,1));
geom.piv.Di = geom.piv.Di(I,:);

clear X I;
