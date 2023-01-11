function fVector = computeF(obj,t,options)
%COMPUTEF Summary of this function goes here
%   Detailed explanation goes here
arguments
    obj ApproxParabolicProblem
    t   
    options.correctConvection logical =true
end

if isnumeric(obj.beta) && all(obj.beta==[0;0])
    options.correctConvection = false ;
end

fVector=zeros(obj.mesh.getDoF(),1);
%for each triangle
for e= 1:obj.mesh.geom.nelements.nTriangles
    %coordinates=obj.mesh.geom.elements.coordinates(obj.mesh.geom.elements.triangles(e,1:3),:);
    %el=Element(coordinates,obj.refElement);
    el=obj.mesh.getEl(obj.refElement,e);
    quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;
    for j=1:el.nDoF
        jj=obj.mesh.geom.pivot.pivot(obj.mesh.geom.elements.triangles(e,j));
        if jj>0

            forcing_term= @(x) obj.f(t,x)* el.phi{j}(x);
            if options.correctConvection
                f_correction=obj.correctForcingConvection(el,j,t);
            else
                f_correction=0;
            end
            fVector(jj) = fVector(jj)...
                +quadrature(forcing_term)+f_correction;
        end

    end


end
end

