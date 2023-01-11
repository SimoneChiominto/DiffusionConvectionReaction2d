function Gn = computeGn(obj,t)
%COMPUTEGN Summary of this function goes here
%   Detailed explanation goes here
Gn=zeros(obj.mesh.getDoF,1);
for e=1:length(obj.mesh.geom.pivot.Ne)
    %take the right boundary function
    g_n=obj.BC.boundaryFunctions{obj.mesh.geom.pivot.Ne(e,2)};

    l=obj.mesh.geom.pivot.Ne(e,1);
    i_b = obj.mesh.geom.elements.borders(l,1);
    i_e= obj.mesh.geom.elements.borders(l,2);
    ii_b=obj.mesh.geom.pivot.pivot(i_b);
    ii_e=obj.mesh.geom.pivot.pivot(i_e);
    if ii_b>0
        Gn(ii_b)= Gn(ii_b)...
            + ((2*g_n(t,[obj.mesh.geom.elements.coordinates(i_b,1);obj.mesh.geom.elements.coordinates(i_b,2)]))...
            +g_n(t,[obj.mesh.geom.elements.coordinates(i_e,1);obj.mesh.geom.elements.coordinates(i_e,2)]))/6 ...
            * norm([obj.mesh.geom.elements.coordinates(i_b,1)-obj.mesh.geom.elements.coordinates(i_e,1), obj.mesh.geom.elements.coordinates(i_b,2)-obj.mesh.geom.elements.coordinates(i_e,2)]);
    end
    if ii_e>0
        Gn(ii_e)= Gn(ii_e)...
            + ((g_n(t,[obj.mesh.geom.elements.coordinates(i_b,1);obj.mesh.geom.elements.coordinates(i_b,2)]))...
            +2*g_n(t,[obj.mesh.geom.elements.coordinates(i_e,1);obj.mesh.geom.elements.coordinates(i_e,2)]))/6 ...
            * norm([obj.mesh.geom.elements.coordinates(i_b,1)-obj.mesh.geom.elements.coordinates(i_e,1), obj.mesh.geom.elements.coordinates(i_b,2)-obj.mesh.geom.elements.coordinates(i_e,2)]);
    end
end %e=1:length(obj.mesh.geom.pivot.Ne)

end

