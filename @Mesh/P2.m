function []=P2(obj)

nnode = obj.geom.nelements.nVertexes;

for l=1:obj.geom.nelements.nBorders

    n(1)=obj.geom.elements.borders(l,1); % edge starting node
    n(2)=obj.geom.elements.borders(l,2); % edge ending node
    e(1)=obj.geom.elements.borders(l,3); % side element
    e(2)=obj.geom.elements.borders(l,4); % side element

    nnode = nnode + 1;

    obj.geom.elements.coordinates(nnode,:) = (obj.geom.elements.coordinates(n(1),:)+...
        obj.geom.elements.coordinates(n(2),:))/2; % coordinates of the edge midpoint

    obj.geom.elements.borders(l,5)=nnode; % to connect the edge with its midpoint

    idx = [1 2 3];

    for el=e

        if(el ~= -1)
            acc = 0;
            acc = idx * ( obj.geom.elements.triangles(el,1:3)==n(1) )';
            acc = acc + idx * ( obj.geom.elements.triangles(el,1:3)==n(2) )';

            switch acc
                case 3
                    obj.geom.elements.triangles(el,4) = nnode;
                case 4
                    obj.geom.elements.triangles(el,6) = nnode;
                case 5
                    obj.geom.elements.triangles(el,5) = nnode;
                otherwise
                    disp('sconoscuto');
            end % switch acc
        end

    end % for el=e

    VertexValue = [0 0];
    Vertex = [0 0];
    D = [0 0];

    %idxV = 1:length(geom.input.BC.InputVertexValues); %indice dei vertici

    % Se il lato e` di bordo
    if( any(e==-1) )
        % Lato di Dirichlet

        if( obj.geom.pivot.nodelist(n(1))~=0 && obj.geom.pivot.nodelist(n(2))~=0 )
            %-----------------------
            if( obj.geom.pivot.nodelist(n(1)) ~= obj.geom.pivot.nodelist(n(2)) )
                %
                if( any(obj.geom.input.BC.InputVertexValues==obj.geom.pivot.nodelist(n(1))) )
                    % Vertice(1) con marker speciale
                    VertexValue(1) = 1;
                    % e` il vettore con un 1 in corrispondenza del vertice
                    % del lato con marker speciale

                    Vertex(1) = obj.geom.input.BC.Boundary.Values*...
                        (obj.geom.input.BC.InputVertexValues == obj.geom.pivot.nodelist(n(1)))';
                    % valore del marker del lato del poligono iniziale che segue n1

                    InputVertexValue = obj.geom.input.BC.InputVertexValues*...
                        (obj.geom.input.BC.InputVertexValues == obj.geom.pivot.nodelist(n(1)))';
                    % indice del nodo del poligono iniziale ce corrisponde al nodo n(1) del mio lato

                    % valore del marker del vertice del poligono iniziale n2
                    D(1) = obj.geom.pivot.nodelist(n(2));
                end
                %
                if( any(obj.geom.input.BC.InputVertexValues==obj.geom.pivot.nodelist(n(2))) )
                    VertexValue(2) = 1;
                    Vertex(2) = obj.geom.input.BC.Boundary.Values*...
                        (obj.geom.input.BC.InputVertexValues == obj.geom.pivot.nodelist(n(2)))';
                    InputVertexValue = obj.geom.input.BC.InputVertexValues*...
                        (obj.geom.input.BC.InputVertexValues == obj.geom.pivot.nodelist(n(2)))';
                    D(2) = obj.geom.pivot.nodelist(n(1));
                end
                %
                if( sum(VertexValue) ~= 2 )
                    Di = VertexValue*D';
                    % nodo con condizione di Dirichlet
                    obj.geom.pivot.nodelist(nnode)= Di;
                    obj.geom.pivot.Di(end+1,:) = [nnode, Di];
                    obj.geom.pivot.pivot(nnode) = min(obj.geom.pivot.pivot)-1;
                else
                    % diamo al nuovo nodo il marker del lato
                    % il lato che stiamo analizzando e` un lato del poligono
                    % iniziale:

                    % l'indice del lato e` quello del nodo di inizio di quel lato
                    if( max(Vertex)-min(Vertex)>1 ) % siamo sul lato di chiusura
                        Di = obj.geom.input.BC.Boundary.Values(max(Vertex));
                    else % siamo sui lati 1->2->3->4->
                        Di = obj.geom.input.BC.Boundary.Values(min(Vertex));
                    end
                    % check della condizione di Neumann aperta
                    if( rem(Di,2)== 0 ) % nodo con grado di liberta`, lato di
                        % Dirichlet aperto
                        obj.geom.pivot.nodelist(nnode)= 0;
                        obj.geom.pivot.pivot(nnode) = max(obj.geom.pivot.pivot)+1;
                    else
                        obj.geom.pivot.nodelist(nnode)= Di;
                        obj.geom.pivot.Di(end+1,:) = [nnode, Di];
                        obj.geom.pivot.pivot(nnode) = min(obj.geom.pivot.pivot)-1;
                    end
                end % if( sum(VertexValue) ~= 2 )
                %----------------------------------
            else % if( geom.pivot.nodelist(n(1)) ~= geom.pivot.nodelist(n(2)) )
                Di = obj.geom.pivot.nodelist(n(1));
                obj.geom.pivot.nodelist(nnode)= Di;
                obj.geom.pivot.Di(end+1,:) = [nnode, Di];
                obj.geom.pivot.pivot(nnode) = min(obj.geom.pivot.pivot)-1;
            end % if( geom.pivot.nodelist(n(1)) ~= geom.pivot.nodelist(n(2)) )
            %----------------------------------

        else
            % Lato di Neumann
            obj.geom.pivot.nodelist(nnode) = 0;
            obj.geom.pivot.pivot(nnode) = max(obj.geom.pivot.pivot)+1;
        end % if( geom.pivot.nodelist(n(1))~=0 & geom.pivot.nodelist(n(2))~=0 )

    else % if( any(e==-1) )
        obj.geom.pivot.nodelist(nnode) = 0;
        obj.geom.pivot.pivot(nnode) = max(obj.geom.pivot.pivot)+1;
    end %if( any(e==-1) )


end % for l=1:geom.nelements.nBorders
end