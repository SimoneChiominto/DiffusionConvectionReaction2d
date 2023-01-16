classdef ApproxDiffusionConvectionReactionProblem2D < DiffusionConvectionReactionProblem2D
    %APPROXDIFFUSIONCONVECTIONREACTIONPROBLEM2D is an handle class that
    % models a the finite element approximation of a Diffusion Convection
    % Reaction Problem in 2D. It is a subclass of DiffusionConvectionReactionProblem2D
    %
    % See also DIFFUSIONCONVECTIONREACTIONPROBLEM2D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties:
    %   - mesh: handle class containing the result of the mesher
    %   - refElement: struct containing the description of the reference
    %     Element
    %   - appproxSolution: is the vector of the approximate value
    %     of the solution computed on the vertices of the mesh
    %   - stiffnessMatrix: is the matrix computed to solve the approximate
    %     problem
    %   - forcingVector: is the vector of known terms for for the solution
    %     of the approximate problem
    %   - error: struct containing the L2 H0 and Linf seminorm


    properties
        mesh
        refElement
        approxSolution
        stiffnessMatrix
        forcingVector
        error
    end

    methods
        function obj = ApproxDiffusionConvectionReactionProblem2D(problem,mesh,refElement)
            %APPROXDIFFUSIONCONVECTIONREACTIONPROBLEM2D Construct an instance of this class
            %   See also APPROXDIFFUSIONCONVECTIONREACTIONPROBLEM2D

            obj=obj@DiffusionConvectionReactionProblem2D(problem.mu,problem.beta,problem.gamma,problem.f,...
                problem.domain)
            obj.BC=problem.BC;
            obj.exactSolution=problem.exactSolution;
            obj.mesh=mesh;
            obj.refElement=refElement;
                if length(obj.refElement.phi)==6
                    obj.mesh.P2()
                end
            obj.error.L2norm=[];
            obj.error.H0seminorm=[];
            obj.error.LInfnorm=[];
        end

        [] = generateLinearSystem(obj,options)
        [diff,conv] = correctStiffnessConvection(obj,el,k,j)
        f = correctForcingConvection(obj,el,j)
        

        function [] = solve(obj)
            %SOLVE solve the linear system corresponding to the degrees of
            %freedom and compute a vector of the approximate solution
            %corresponding to all vertices of the mesh

            if isempty(obj.stiffnessMatrix)
                obj.generateLinearSystem()
            end

            %if stiffness matrix is symmetric we can use pcg
            if all(obj.stiffnessMatrix==obj.stiffnessMatrix')
                sol = pcg(obj.stiffnessMatrix,obj.forcingVector,1.0E-8,500);
            else
                sol=obj.stiffnessMatrix\obj.forcingVector;
            end
            
            %add solution to vertices that are not degrees of freedom
            
            obj.approxSolution=zeros(length(obj.mesh.geom.elements.coordinates),1);
            for j=1:length(obj.mesh.geom.elements.coordinates)
                jj=obj.mesh.geom.pivot.pivot(j);
                if jj>0
                    obj.approxSolution(j)=sol(jj);
                else
                    g_d=obj.BC.boundaryFunctions{obj.mesh.geom.pivot.Di(-jj,2)};
                    obj.approxSolution(j)=g_d(obj.mesh.geom.elements.coordinates(j,:)');
                end
            end

        end



        function norm= getL2Error(obj)
            %GETL2ERROR compute L2 error of the approximate solution
            if isempty(obj.exactSolution.u)
                error("no exact solution known")
            end
            if isempty(obj.error.L2norm)
                norm=0;
                %quadrature_ref= @(f) 1/6*(f([0;1/2])+f([1/2;0])+f([1/2;1/2]));

                for e= 1:obj.mesh.geom.nelements.nTriangles
                    coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,:),:);
                    el=Element(coordinates,obj.refElement);
                    quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

                    Ge=obj.mesh.geom.elements.triangles(e,:);

                    u_approx= @(x) dot(obj.approxSolution(Ge),EvalCell(el.phi,x));
                    %obj.approxSolution(Ge(1))*el.phi{1}(x) ...
                    %    + obj.approxSolution(Ge(2))*el.phi{2}(x)...
                    %    + obj.approxSolution(Ge(3))*el.phi{3}(x);
                    norm=norm+quadrature (@(x) (obj.exactSolution.u(x)-u_approx(x))^2);

                end

                norm=norm.^0.5;
                obj.error.L2norm=norm;
                return
            end

            norm=obj.error.L2norm;
        end

        function norm=getH0Error(obj)
            %GETH0ERROR compute H0 semierror of the approximate solution
            if isempty(obj.exactSolution.uGrad)
                error("no gradient of exact solution known")
            end
            if isempty(obj.error.H0seminorm)
                norm=0;
                %quadrature_ref= @(f) 1/6*(f([0;1/2])+f([1/2;0])+f([1/2;1/2]));

                for e= 1:obj.mesh.geom.nelements.nTriangles

                    coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,:),:);
                    el=Element(coordinates,obj.refElement);
                    quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

                    Ge=obj.mesh.geom.elements.triangles(e,:);

                    grad_u_approx= @(x) EvalCell(el.gradPhi,x)*obj.approxSolution(Ge);
                    %el.gradPhi{1}(x) ...
                    %    + obj.approxSolution(Ge(2))*el.gradPhi{2}(x) ...
                    %    + obj.approxSolution(Ge(3))*el.gradPhi{3}(x);

                    norm=norm+quadrature(@(x) sum((obj.exactSolution.uGrad(x) - grad_u_approx(x)).^2) );
                end

                norm=norm^0.5;
                obj.error.H0seminorm=norm;
                return
            end

            norm=obj.error.H0seminorm;

        end

        function norm=getLInfError(obj)
            %GETLINFERROR compute Linf error of the approximate solution
            if isempty(obj.exactSolution.u)
                error("no exact solution known")
            end
            if isempty(obj.error.LInfnorm)
                norm=max(abs(obj.approxSolution-obj.exactSolution.u(obj.mesh.geom.elements.coordinates')'));
                obj.error.LInfnorm=norm;
                return
            end

            norm=obj.error.LInfnorm;
        end



        function fig=plot(obj,options)
            arguments
                obj ApproxDiffusionConvectionReactionProblem2D
                options.type ="solution"
            end
            if isempty(obj.approxSolution)
                obj.solve();
            end
            fig=figure;

            if options.type=="solution"
                z=obj.approxSolution;
                color_map="autumn";
            elseif options.type=="error"
                if isempty(obj.exactSolution.u)
                    error("No exact solution known")
                end
                z=obj.exactSolution.u(obj.mesh.geom.elements.coordinates')'-obj.approxSolution;
                color_map="winter";
            elseif options.type=="exact_solution"
                if isempty(obj.exactSolution.u)
                    error("No exact solution known")
                end
                z=obj.exactSolution.u(obj.mesh.geom.elements.coordinates')';
                color_map="summer";
            else
                error("type not known, try asking solution or error")
            end

            pdeplot(obj.mesh.geom.elements.coordinates',...
                obj.mesh.geom.elements.triangles',...
                "ZData",z,...
                "XYData",z, ...
                "Mesh","on","ColorMap",color_map,...
                "XYStyle","flat")
        end

    end
end


function out=EvalCell(cellArray,x)
    out=zeros(length(cellArray{1}(x)),length(cellArray));
    for i=1:length(cellArray)
        out(:,i)=cellArray{i}(x);
        
    end
end

