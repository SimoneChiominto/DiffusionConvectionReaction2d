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

%quadrature_ref= @(f) 1/6*(f([0;1/2])+f([1/2;0])+f([1/2;1/2]));

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
                    [diff_correction,conv_correction]=obj.correctStiffnessConvection(el,k,j);
                    obj.stiffnessMatrix(jj,kk)= obj.stiffnessMatrix(jj,kk)  ...
                        +quadrature(diffusion)+quadrature(convection)...%+quadrature(reaction)...
                        +diff_correction+conv_correction;
                    obj.stiffnessMatrix(jj,jj)=obj.stiffnessMatrix(jj,jj)+quadrature(reaction);

                else
                    [diff_correction,conv_correction]=obj.correctStiffnessConvection(el,k,j);
                    Ad(jj,-kk)= Ad(jj,-kk) +...
                        +quadrature(diffusion)+quadrature(convection)...%+quadrature(reaction)...
                        +diff_correction+conv_correction;
                    obj.stiffnessMatrix(jj,jj)=obj.stiffnessMatrix(jj,jj)+quadrature(reaction);
                    %countD=countD+1;
                    %JJ(countD)=jj;
                    %KK(countD)=-kk;
                    %val(count)=quadrature(diffusion)+quadrature(convection)+quadrature(reaction);

                end %if kk>0
            end %for k=1:3

            forcing_term= @(x) f(x)* el.phi{j}(x);
            f_correction=obj.correctForcingConvection(el,j);
            b(jj) = b(jj)...
                +quadrature(forcing_term)+f_correction;


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

