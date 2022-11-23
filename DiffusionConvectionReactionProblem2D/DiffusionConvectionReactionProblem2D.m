classdef DiffusionConvectionReactionProblem2D < handle
    % DiffusionConvectionReactionProblem2D is an handle class that models a 
    % Diffusion Convection Reaction Problem in 2D 
    %   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Properties:
    %   - mu: diffusion function handle
    %   - beta: convection function handle
    %   - gamma: reaction function handle
    %   - f: forcing term function handle
    %   - domain: is a struct containg a description of the domain where
    %     the problem is defined
    %   - BC: is a struct containing a description of the boundary
    %     conditions of the problem
    %   - exactSolution: is a struct with two attributes:
    %       - u: the exact solution function handle
    %       - uGrad: function handle containing the gradient of the exact solution
    
    properties
        mu              %diffusion
        beta            %convection
        gamma           %reaction
        f               %forcing term
        domain          %Domain defined as in the Mesher input
        BC              %BC defined as in the Mesher input with some other properties
        exactSolution   %exact solution and gradient of the exact solution
    end
    
    methods
        function obj = DiffusionConvectionReactionProblem2D(mu,beta,gamma,f,domain,BC,boundaryFunctions, u, uGrad)
            %DiffusionConvectionReactionProblem2D: Construct an instance
            %the class DiffusionConvectionReactionProblem2D
            %
            %See also DIFFUSIONCONVECTIONREACTIONPROBLEM2D  
            if nargin==0
                return
            end
            obj.mu = mu;
            obj.beta = beta;
            obj.gamma = gamma;
            obj.f = f;
            obj.domain = domain;
            if nargin>=7
                obj.BC=DefineBoundaryConditions(obj,BC,boundaryFunctions);
            end
            if nargin>=9
                obj.exactSolution.u = u;
                obj.exactSolution.uGrad = uGrad;                
            end            
        end %DiffusionConvectionReactionProblem2D
        
        function BC = DefineBoundaryConditions(obj,BC,boundaryFunctions)
            %DEFINEBOUNDARYCONDITIONS create a cell array of function 
            % handles, which represent the functions associated to a boundary 
            %   
            %   This function modifies the cell array boundaryFunctions adding
            %   cells corresponding to vertices of the domain the corresponds
            %   to dirichlet conditions.
            %   If a vertices is in the intersection of two sides
            %   corresponding to dirichlet condition this function check if
            %   the function is continuous.
            %   It add this modyfied boundaryFunctions to the BC struct
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Inputs: 
            %
            %   BC: BC is a struct containing informations on the Boundary
            %       Conditions, 
            %       See also DIFFUSIONCONVECTIONREACTIONPROBLEM2D
            %
            %   boundaryFunctions: boundaryFunctions is a cell array
            %       boundaryFunctions{k} contains the function that
            %       corresponds to the side of the polygon associated to the
            %       marker j;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            N_vert=length(BC.InputVertexValues);
            %add a new component to the cell array for each verices
            for v=1:N_vert
                % see if the sides of the vertices are dirichlet or neumann
                % sides
                vert_conds=[BC.Boundary.Values(obj.domain.Boundary.Values(mod(v-2,N_vert)+1)),BC.Boundary.Values(obj.domain.Boundary.Values(mod(v-1,N_vert)+1))];
                if all(mod(vert_conds,2)) 
                    %check if function is continuous
                    if boundaryFunctions{vert_conds(1)}(obj.domain.InputVertex(v,:)')~=boundaryFunctions{vert_conds(2)}(obj.domain.InputVertex(v,:)')
                        error("Function on the border is not continuous")
                    end
                end
                for edge=1:2
                    if mod(vert_conds(edge),2)
                        boundaryFunctions{BC.InputVertexValues(v)}=boundaryFunctions{vert_conds(edge)};
                        break
                    end
                end

            end
            
            BC.boundaryFunctions = boundaryFunctions;
            
        end %DefineBoundaryConditions

        function fig=plotSolution(obj,mesh)
            %PLOTSOLUTION plots the exact solution of a problem

            if isempty(obj.exactSolution)
                error("no exact solution known")
            end
            fig=figure;
            
            if nargin==1 
                fsurf(@(x,y) obj.exactSolution.u([x;y]), ...
                    [min(obj.domain.InputVertex(:,1)),...
                    max(obj.domain.InputVertex(:,1)),...
                    min(obj.domain.InputVertex(:,2)),...
                    max(obj.domain.InputVertex(:,2))]);
                return
            end

            trisurf(mesh.geom.elements.triangles, ...
                mesh.geom.elements.coordinates(:,1),...
                mesh.geom.elements.coordinates(:,2),...
                obj.exactSolution.u(mesh.geom.elements.coordinates'));
        end %plotSolution


    end
end

