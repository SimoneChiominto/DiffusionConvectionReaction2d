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

        function [] = generateLinearSystem(obj)
            %GENERATELINEARSYSTEM Compute Stiffness matrix and forcing
            %vector

            %save useful variables
            f=obj.f;
            mu=obj.mu;
            beta=obj.beta;
            gamma=obj.gamma;
            boundaryFunctions=obj.BC.boundaryFunctions;
            n_dof=obj.mesh.getDoF();
            n_d=-min(obj.mesh.geom.pivot.pivot);
            
            %preallocate vectors ad matrices
            b=zeros(n_dof,1);
            b_n=zeros(n_dof,1);
            obj.stiffnessMatrix=spalloc(n_dof,n_dof,10*n_dof);
            Ad=spalloc(n_dof,n_d,10*n_dof);
            %JJ=zeros(10*n_dof,1);
            %KK=zeros(10*n_dof,1);
            %val=zeros(10*n_dof,1);
            %JJD=zeros(10*n_dof,1);
            %KKD=zeros(10*n_dof,1);
            %valD=zeros(10*n_dof,1);
            %count=0;
            %countD=0;

            quadrature_ref= @(f) 1/6*(f([0;1/2])+f([1/2;0])+f([1/2;1/2]));

            %for each triangle
            for e= 1:obj.mesh.geom.nelements.nTriangles
                
                %create the real element object
                coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,1:3),:);
                el=Element(coordinates,obj.refElement);
                quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

                %for each vertex of the triangle
                for j=1:el.nDoF
                    jj=obj.mesh.geom.pivot.pivot(obj.mesh.geom.elements.triangles(e,j));
                    
                    %check if it is a degree of freedom
                    if jj>0 
                            
                        %for each vertex of the triangle
                        for k=1:el.nDoF
                            kk=obj.mesh.geom.pivot.pivot(obj.mesh.geom.elements.triangles(e,k));
                            
                            diffusion= @(x) mu(x) .* el.gradPhi{k}(x)'*el.gradPhi{j}(x);
                            convection = @(x) beta(x)' * el.gradPhi{k}(x) * el.phi{j}(x);
                            reaction= @(x) gamma(x) * el.phi{k}(x) * el.phi{j}(x);

                            %check if it is a degree of freedom
                            if kk>0
                                %count=count+1;
                                %JJ(count)=jj;
                                %KK(count)=kk;
                                %val(count)=quadrature(diffusion)+quadrature(convection)+quadrature(reaction);
                                obj.stiffnessMatrix(jj,kk)= obj.stiffnessMatrix(jj,kk)  ...
                                    +quadrature(diffusion)+quadrature(convection)+quadrature(reaction);

                            else
                                Ad(jj,-kk)= Ad(jj,-kk) +...
                                    +quadrature(diffusion)+quadrature(convection)+quadrature(reaction);
                                %countD=countD+1;
                                %JJ(countD)=jj;
                                %KK(countD)=-kk;
                                %val(count)=quadrature(diffusion)+quadrature(convection)+quadrature(reaction);

                            end %if kk>0
                        end %for k=1:3

                        forcing_term= @(x) f(x)* el.phi{j}(x);
                        b(jj) = b(jj)...
                            +quadrature(forcing_term);


                    end %if jj>0

                end %for j=1:3
            end %for e= 1:obj.mesh.geom.nelements.nTriangles

            
            for e=1:length(obj.mesh.geom.pivot.Ne)
                %take the right boundary function
                g_n=boundaryFunctions{obj.mesh.geom.pivot.Ne(e,2)};
                
                l=obj.mesh.geom.pivot.Ne(e);
                i_b = obj.mesh.geom.elements.borders(l,1);
                i_e= obj.mesh.geom.elements.borders(l,2);
                ii_b=obj.mesh.geom.pivot.pivot(i_b);
                ii_e=obj.mesh.geom.pivot.pivot(i_e);
                if ii_b>0
                    b_n(ii_b)= b_n(ii_b)...
                        + ((2*g_n([obj.mesh.geom.elements.coordinates(i_b,1);obj.mesh.geom.elements.coordinates(i_b,2)]))+g_n([obj.mesh.geom.elements.coordinates(i_e,1);obj.mesh.geom.elements.coordinates(i_e,2)]))/6 ...
                        * norm([obj.mesh.geom.elements.coordinates(i_b,1)-obj.mesh.geom.elements.coordinates(i_e,1), obj.mesh.geom.elements.coordinates(i_b,2)-obj.mesh.geom.elements.coordinates(i_e,2)]);
                end
                if ii_e>0
                    b_n(ii_e)= b_n(ii_e)...
                        + ((g_n([obj.mesh.geom.elements.coordinates(i_b,1);obj.mesh.geom.elements.coordinates(i_b,2)]))+2*g_n([obj.mesh.geom.elements.coordinates(i_e,1);obj.mesh.geom.elements.coordinates(i_e,2)]))/6 ...
                        * norm([obj.mesh.geom.elements.coordinates(i_b,1)-obj.mesh.geom.elements.coordinates(i_e,1), obj.mesh.geom.elements.coordinates(i_b,2)-obj.mesh.geom.elements.coordinates(i_e,2)]);
                end
            end %e=1:length(obj.mesh.geom.pivot.Ne)

            u_d=zeros(length(obj.mesh.geom.pivot.Di),1);
            for i=1:length(obj.mesh.geom.pivot.Di)
                %take the right boundary function
                g_d=boundaryFunctions{obj.mesh.geom.pivot.Di(i,2)};
                u_d(i)=g_d(obj.mesh.geom.elements.coordinates(obj.mesh.geom.pivot.Di(i,1),:)');
            end
            
            %compute forcing vector
            obj.forcingVector=b-Ad*u_d+b_n;

        end %generateLinearSystem


        function [] = solve(obj)
            %SOLVE solve the linear system corresponding to the degrees of
            %freedom and compute a vector of the approximate solution
            %corresponding to all vertices of the mesh

            if isempty(obj.stiffnessMatrix)
                obj.generateLinearSystem()
            end

            %if stiffness matrix is symmetric we can use pcg
            if all(obj.stiffnessMatrix==obj.stiffnessMatrix')
                sol = pcg(obj.stiffnessMatrix,obj.forcingVector,0.001,200);
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
                quadrature_ref= @(f) 1/6*(f([0;1/2])+f([1/2;0])+f([1/2;1/2]));

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
                quadrature_ref= @(f) 1/6*(f([0;1/2])+f([1/2;0])+f([1/2;1/2]));

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
                norm=max(abs(obj.approxSolution-obj.exactSolution.u(obj.mesh.geom.elements.coordinates')));
                obj.error.LInfnorm=norm;
                return
            end

            norm=obj.error.LInfnorm;
        end



        function fig=plot(obj)
            if isempty(obj.approxSolution)
                obj.solve();
            end
            fig=figure;
            trisurf(obj.mesh.geom.elements.triangles,...
                obj.mesh.geom.elements.coordinates(:,1),...
                obj.mesh.geom.elements.coordinates(:,2),...
                obj.approxSolution)
        end


    end
end


function out=EvalCell(cellArray,x)
    out=zeros(length(cellArray{1}(x)),length(cellArray));
    for i=1:length(cellArray)
        out(:,i)=cellArray{i}(x);
    end
end

