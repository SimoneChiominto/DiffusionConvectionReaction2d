function [] = computeMassMatrix(obj,options)
%COMPUTEMASSMATRIX Summary of this function goes here
%   Detailed explanation goes here
arguments
    obj ApproxParabolicProblem
    options.massLumping         logical = false
end


n_dof=obj.mesh.getDoF();
n_d=-min(obj.mesh.geom.pivot.pivot);
JJ=zeros((2*options.massLumping+1)*18*n_dof,1);
KK=zeros((2*options.massLumping+1)*18*n_dof,1);
val=zeros((2*options.massLumping+1)*18*n_dof,1);
count=1;
JJ_D=zeros(10*n_dof,1);
KK_D=zeros(10*n_dof,1);
val_D=zeros(10*n_dof,1);
count_D=1;

for e= 1:obj.mesh.geom.nelements.nTriangles

    %coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,1:3),:);
    %el=Element(coordinates,obj.refElement);
    el=obj.mesh.getEl(obj.refElement,e);
    
    quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;
    for j=1:el.nDoF
        jj=obj.mesh.geom.pivot.pivot(obj.mesh.geom.elements.triangles(e,j));

        %check if it is a degree of freedom
        if jj>0
        
            %for each vertex of the triangle
            for k=1:el.nDoF
                kk=obj.mesh.geom.pivot.pivot(obj.mesh.geom.elements.triangles(e,k));        
                
                if el.type~="P1"
                    mass= @(x) el.phi{k}(x) * el.phi{j}(x);
                end

                if kk>0
                    if el.type=="P1"
                        mass_component=(el.Area*(1+(j==k)))/12;
                    else
                        mass_component=quadrature(mass);
                    end

                    JJ(count)=jj;
                    KK(count)=kk;
                    val(count)=(~options.massLumping)* mass_component;
                    count=count+(val(count)~=0)*1;
                    
                    if options.massLumping
                        JJ(count)=jj;
                        KK(count)=jj;
                        val(count)=mass_component;
                        count=count+(val(count)~=0)*1;
                    end
                  
                else
                    if el.type=="P1"
                        mass_component=(el.Area*(1+(j==k)))/12;
                    else
                        mass_component=quadrature(mass);
                    end

                    JJ_D(count_D)=jj;
                    KK_D(count_D)=-kk;
                    val_D(count_D)=(~options.massLumping)* mass_component;
                    count_D=count_D+(val_D(count_D)~=0)*1;

                    if options.massLumping
                        JJ(count)=jj;
                        KK(count)=jj;
                        val(count)=reaction_component;
                        count=count+(val(count)~=0)*1;
                    end
                    
                end
            end %k=1:el.nDoF
        end %if jj>0

    end %j=1:el.nDoF
end %for e= 1:obj.mesh.geom.nelements.nTriangles
obj.massMatrix=sparse(JJ(1:count-1),KK(1:count-1),val(1:count-1),n_dof,n_dof);
obj.Md=sparse(JJ_D(1:count_D-1),KK_D(1:count_D-1),val_D(1:count_D-1),n_dof,n_d);

end

