function [] = computeStiffnessMatrix(obj,options)
%COMPUTESTIFFNESSMATRIX Summary of this function goes here
%   Detailed explanation goes here
arguments
    obj ApproxParabolicProblem
    options.massLumping         logical = false
    options.correctConvection   logical = true
end

if isnumeric(obj.beta) && all(obj.beta==[0;0])
    options.correctConvection = false ;
end

n_dof=obj.mesh.getDoF();
n_d=-min(obj.mesh.geom.pivot.pivot);

JJ=zeros((2*options.massLumping+1)*18*n_dof,1);
KK=zeros((2*options.massLumping+1)*18*n_dof,1);
val=zeros((2*options.massLumping+1)*18*n_dof,1);
JJ_D=zeros(10*n_dof,1);
KK_D=zeros(10*n_dof,1);
val_D=zeros(10*n_dof,1);
count=1;
countD=1;

%for each triangle
for e= 1:obj.mesh.geom.nelements.nTriangles

    %create the real element object
    %coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,1:3),:);
    %el=Element(coordinates,obj.refElement);
    el=obj.mesh.getEl(obj.refElement,e);
    
    if any([isnumeric(obj.mu),isnumeric(obj.beta),isnumeric(obj.gamma)]) && el.type=="P1"
        dx=el.getdx();
        dy=el.getdy();
    end

    quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

    %for each vertex of the triangle
    for j=1:el.nDoF
        jj=obj.mesh.geom.pivot.pivot(obj.mesh.geom.elements.triangles(e,j));

        %check if it is a degree of freedom
        if jj>0

            %for each vertex of the triangle
            for k=1:el.nDoF
                kk=obj.mesh.geom.pivot.pivot(obj.mesh.geom.elements.triangles(e,k));        
                    
                if isnumeric(obj.mu) && el.type~="P1"
                    diffusion= @(x) obj.mu .* el.gradPhi{k}(x)'*el.gradPhi{j}(x);
                elseif isa(obj.mu,"function_handle")
                    diffusion= @(x) obj.mu(x) .* el.gradPhi{k}(x)'*el.gradPhi{j}(x);
                end

                if isnumeric(obj.beta) && el.type~="P1"
                    convection = @(x) obj.beta' * el.gradPhi{k}(x) * el.phi{j}(x);
                elseif isa(obj.beta,"function_handle")
                    convection = @(x) obj.beta(x)' * el.gradPhi{k}(x) * el.phi{j}(x);
                end

                if isnumeric(obj.gamma) && el.type~="P1"
                    reaction= @(x) obj.gamma * el.phi{k}(x) * el.phi{j}(x);
                elseif isa(obj.gamma,"function_handle")
                    reaction= @(x) obj.gamma(x) * el.phi{k}(x) * el.phi{j}(x);
                end
                %check if it is a degree of freedom
                
                if kk>0
                    if options.correctConvection
                        [diff_correction,conv_correction]=obj.correctStiffnessConvection(el,k,j);
                    else
                        diff_correction=0;
                        conv_correction=0;
                    end

                    if isnumeric(obj.mu) && el.type=="P1"
                        diffusion_component= obj.mu * (dy(k)*dy(j)+dx(k)*dx(j))/(4*el.Area);
                    else
                        diffusion_component=quadrature(diffusion);
                    end

                    if isnumeric(obj.beta) && el.type=="P1"
                        convection_component= (obj.beta(1)*dy(k)+obj.beta(2)*dx(k))/6;
                    else
                        convection_component=quadrature(convection);
                    end

                    if isnumeric(obj.gamma) && el.type=="P1"
                        reaction_component=(obj.gamma*el.Area*(1+(j==k)))/12;
                    else
                        reaction_component=quadrature(reaction);
                    end
                    
                    
                    JJ(count)=jj;
                    KK(count)=kk;
                    val(count)=diffusion_component+convection_component...
                        +(~options.massLumping)* reaction_component...
                        +diff_correction+conv_correction;
                    count=count+(val(count)~=0)*1;
                    
                    if options.massLumping
                        JJ(count)=jj;
                        KK(count)=jj;
                        val(count)=reaction_component;
                        count=count+(val(count)~=0)*1;
                    end
                    %obj.stiffnessMatrix(jj,jj)=obj.stiffnessMatrix(jj,jj)+quadrature(reaction);

                else
                    if options.correctConvection
                        [diff_correction,conv_correction]=obj.correctStiffnessConvection(el,k,j);
                    else
                        diff_correction=0;
                        conv_correction=0;
                    end

                    if isnumeric(obj.mu) && el.type=="P1"
                        diffusion_component= obj.mu * (dy(k)*dy(j)+dx(k)*dx(j))/(4*el.Area);
                    else
                        diffusion_component=quadrature(diffusion);
                    end

                    if isnumeric(obj.beta) && el.type=="P1"
                        convection_component= (obj.beta(1)*dy(k)+obj.beta(2)*dx(k))/6;
                    else
                        convection_component=quadrature(convection);
                    end

                    if isnumeric(obj.gamma) && el.type=="P1"
                        reaction_component=(obj.gamma*el.Area*(1+(j==k)))/12;
                    else
                        reaction_component=quadrature(reaction);
                    end
                    
                    
                    JJ_D(countD)=jj;
                    KK_D(countD)=-kk;
                    val_D(countD)=diffusion_component+convection_component...
                        +(~options.massLumping)* reaction_component...
                        +diff_correction+conv_correction;
                    countD=countD+(val_D(countD)~=0)*1;
                    
                    if options.massLumping
                        JJ(count)=jj;
                        KK(count)=jj;
                        val(count)=reaction_component;
                        count=count+(val(count)~=0)*1;
                    end
                end %if kk>0
            
            end %for k=1:el.nDoF
        
        end %if jj>0

    end %for j=1:el.nDoF

end %for e= 1:obj.mesh.geom.nelements.nTriangles

obj.stiffnessMatrix=sparse(JJ(1:count-1),KK(1:count-1),val(1:count-1),n_dof,n_dof);
obj.Ad=sparse(JJ_D(1:countD-1),KK_D(1:countD-1),val_D(1:countD-1),n_dof,n_d);

end

