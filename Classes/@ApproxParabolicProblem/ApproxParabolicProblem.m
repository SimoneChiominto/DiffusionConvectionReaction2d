classdef ApproxParabolicProblem < ParabolicProblem
    %APPROXPARABOLICPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        timeDiscretization
        refElement
        approxSolution
        stiffnessMatrix
        Ad
        massMatrix
        Md
        forcingVector
        Gn
        error
    end

    
    methods
        function obj = ApproxParabolicProblem(problem,mesh,refElement,timeDiscretization)
            %APPROXDIFFUSIONCONVECTIONREACTIONPROBLEM2D Construct an instance of this class
            %   See also APPROXDIFFUSIONCONVECTIONREACTIONPROBLEM2D

            obj=obj@ParabolicProblem(problem.mu,problem.beta,problem.gamma,problem.f,...
                problem.domain,problem.timeInterval);
            obj.BC=problem.BC;
            obj.exactSolution=problem.exactSolution;
            obj.mesh=mesh;
            if length(timeDiscretization)>1
                obj.timeDiscretization=timeDiscretization;
            else
                obj.timeDiscretization=linspace(problem.timeInterval(1),problem.timeInterval(2),timeDiscretization+1);
            end

            obj.refElement=refElement;
                if length(obj.refElement.phi)==6
                    obj.mesh.P2()
                end
            obj.error.L2norm=[];
            obj.error.H0seminorm=[];
            obj.error.LInfnorm=[];
        end
 
        [] = computeMassMatrix(obj,options)
        [] = computeStiffnessMatrix(obj,options)
        [diff,conv] = correctStiffnessConvection(obj,el,k,j)
        Gn = computeGn(obj,t)
        fVector = computeF(obj,t)
        f = correctForcingConvection(obj,el,j,n)


        function []=solve(obj)
            arguments
                obj ApproxParabolicProblem
            end

            obj.approxSolution=zeros(length(obj.mesh.geom.elements.coordinates),length(obj.timeDiscretization));
            u0=zeros(obj.mesh.getDoF(),length(obj.timeDiscretization));
            uD=zeros(length(obj.mesh.geom.pivot.Di),length(obj.timeDiscretization));
            
            obj.approxSolution(:,1)=obj.BC.initialSolution(obj.mesh.geom.elements.coordinates')';
            
            
            for j=1:length(obj.mesh.geom.elements.coordinates)
                jj=obj.mesh.geom.pivot.pivot(j);
                if jj>0           
                    u0(jj,1)=obj.BC.initialSolution(obj.mesh.geom.elements.coordinates(j,:)');
                end
            end

            for i=1:length(obj.mesh.geom.pivot.Di)
                uD(i,1)=obj.BC.initialSolution(obj.mesh.geom.elements.coordinates(obj.mesh.geom.pivot.Di(i,1),:)');
            end
            fig=waitbar(0);
            obj.Gn(:,1)=obj.computeGn(obj.timeDiscretization(1));
            obj.forcingVector(:,1)=obj.computeF(obj.timeDiscretization(1));
            
            max_time=length(obj.timeDiscretization);
            for n=1:max_time-1
                waitbar(n/max_time,fig,sprintf("%d di %d",n+1,max_time))
                %compute uD al tempo t_{n+1}
                for i=1:length(obj.mesh.geom.pivot.Di)
                    %take the right boundary function
                    g_d=obj.BC.boundaryFunctions{obj.mesh.geom.pivot.Di(i,2)};
                    uD(i,n+1)=g_d(obj.timeDiscretization(n+1),obj.mesh.geom.elements.coordinates(obj.mesh.geom.pivot.Di(i,1),:)');
                end

                dt=obj.timeDiscretization(n+1)-obj.timeDiscretization(n);
                
                obj.Gn(:,n+1)=obj.computeGn(obj.timeDiscretization(n+1));
                obj.forcingVector(:,n+1)=obj.computeF(obj.timeDiscretization(n+1));
                
                A=(obj.massMatrix+dt/2*obj.stiffnessMatrix);
                b=obj.massMatrix*u0(:,n)...
                    - obj.Md*(uD(:,n+1)-uD(:,n))...
                    -dt/2*obj.stiffnessMatrix*u0(:,n)...
                    -dt/2*obj.Ad*(uD(:,n+1)+uD(:,n))...
                    +dt/2*(obj.forcingVector(:,n)+obj.forcingVector(:,n+1))...
                    +dt/2*(obj.Gn(:,n)+obj.Gn(:,n+1));

                u0(:,n+1)=A\b;
            end
            
            
            for j=1:length(obj.mesh.geom.elements.coordinates)
                jj=obj.mesh.geom.pivot.pivot(j);
                if jj>0           
                    obj.approxSolution(j,:)=u0(jj,:);
                else
                    obj.approxSolution(j,:)=uD(-jj,:);
                end
            end

        end

        function fig=plot(obj,n,fig)
            if isempty(obj.approxSolution)
                obj.solve();
            end
            if nargin<3
                fig=figure;
            end
            %n_vert=obj.mesh.geom.nelements.nVertexes;
            %trisurf(obj.mesh.geom.elements.triangles(:,1:3),...
                %obj.mesh.geom.elements.coordinates(1:n_vert,1),...
                %obj.mesh.geom.elements.coordinates(1:n_vert,2),...
                %obj.approxSolution(1:n_vert))
            if nargin<=1
                pdeplot(obj.mesh.geom.elements.coordinates',...
                    obj.mesh.geom.elements.triangles',...
                    "ZData",obj.approxSolution(:,end),...
                    "XYData",obj.approxSolution(:,end), ...
                    "Mesh","on","ColorMap","autumn",...
                    "XYStyle","flat")
                zlim([min(obj.approxSolution,[],'all'),max(obj.approxSolution,[],'all')]);
            else
                
                pdeplot(obj.mesh.geom.elements.coordinates',...
                    obj.mesh.geom.elements.triangles',...
                    "ZData",obj.approxSolution(:,n),...
                    "XYData",obj.approxSolution(:,n), ...
                    "Mesh","on","ColorMap","autumn",...
                    "XYStyle","flat")
                zlim([min(obj.approxSolution,[],'all'),max(obj.approxSolution,[],'all')]);
            end
        end

        function M=video(obj)
            fig=figure;
            fig.Visible = "off";
            for n=1:length(obj.timeDiscretization)
                plot(obj,n,fig);
                M(n)=getframe;
            end
            fig.Visible = "on";
            movie(M)
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
                    %coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,:),:);
                    %el=Element(coordinates,obj.refElement);
                    el=obj.mesh.getEl(obj.refElement,e);
                    quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

                    Ge=obj.mesh.geom.elements.triangles(e,:);

                    u_approx= @(x) dot(obj.approxSolution(Ge',end),EvalCell(el.phi,x));
                    %obj.approxSolution(Ge(1))*el.phi{1}(x) ...
                    %    + obj.approxSolution(Ge(2))*el.phi{2}(x)...
                    %    + obj.approxSolution(Ge(3))*el.phi{3}(x);
                    norm=norm+quadrature (@(x) (obj.exactSolution.u(obj.timeInterval(2),x)-u_approx(x))^2);

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

                    %coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,:),:);
                    %el=Element(coordinates,obj.refElement);
                    el=obj.mesh.getEl(obj.refElement,e);
                    quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

                    Ge=obj.mesh.geom.elements.triangles(e,:);

                    grad_u_approx= @(x) EvalCell(el.gradPhi,x)*obj.approxSolution(Ge',end);
                    %el.gradPhi{1}(x) ...
                    %    + obj.approxSolution(Ge(2))*el.gradPhi{2}(x) ...
                    %    + obj.approxSolution(Ge(3))*el.gradPhi{3}(x);

                    norm=norm+quadrature(@(x) sum((obj.exactSolution.uGrad(obj.timeInterval(2),x) - grad_u_approx(x)).^2) );
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
                norm=max(abs(obj.approxSolution(:,end)-obj.exactSolution.u(obj.timeInterval(2),obj.mesh.geom.elements.coordinates')'));
                obj.error.LInfnorm=norm;
                return
            end

            norm=obj.error.LInfnorm;
        end
        




    end
end

function out=EvalCell(cellArray,x)
    out=zeros(length(cellArray{1}(x)),length(cellArray));
    for i=1:length(cellArray)
        out(:,i)=cellArray{i}(x);
        
    end
end

