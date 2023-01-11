classdef ParabolicProblem < handle
    %PARABOLICPROBLEM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mu      {mustBeA(mu,["numeric","function_handle"])} %diffusion
        beta    {mustBeA(beta,["numeric","function_handle"])} %convection
        gamma   {mustBeA(gamma,["numeric","function_handle"])}       %reaction
        f       {mustBeA(f,["numeric","function_handle"])}        %forcing term
        domain          %Domain defined as in the Mesher input
        timeInterval
        BC              %BC defined as in the Mesher input with some other properties
        exactSolution   %exact solution and gradient of the exact solution
    end


    methods
        function obj = ParabolicProblem(mu,beta,gamma,f,domain,timeInterval,BC,boundaryFunctions,initialSolution, u, uGrad)
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
            obj.timeInterval=timeInterval;
            if nargin>=8
                obj.BC=DefineBoundaryConditions(obj,BC,boundaryFunctions);
                obj.BC.initialSolution=initialSolution;
            end
            if nargin>=11
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
                    if boundaryFunctions{vert_conds(1)}(0,obj.domain.InputVertex(v,:)')~=boundaryFunctions{vert_conds(2)}(0,obj.domain.InputVertex(v,:)')
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

        


        
    end
end

