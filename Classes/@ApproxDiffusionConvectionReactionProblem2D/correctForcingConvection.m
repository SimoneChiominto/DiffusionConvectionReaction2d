function f = correctForcingConvection(obj,el,j)
%CONVECTIONCORRECTION Summary of this function goes here
%   Detailed explanation goes here


if length(obj.mesh.geom.elements.triangles(1,:))==3
    m_k = 1/3;
elseif length(obj.mesh.geom.elements.triangles(1,:))==6
    m_k = 1/24;
end

quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

h_E=el.getMaxLength();

if isnumeric(obj.beta)
    norm_beta_E = norm(obj.beta);
else
    norm_beta_E = quadrature(@(x) norm(obj.beta(x))) / el.Area ;
end

if isnumeric(obj.mu)
    norm_mu_E = obj.mu;
else
    norm_mu_E = quadrature(@(x) obj.mu(x)) / el.Area;
end

Pe_h = m_k*norm_beta_E*h_E / (2*norm_mu_E);

if Pe_h<=1
    tau_E=m_k*(h_E^2)/(4*norm_mu_E);
else
    tau_E=h_E/(2*norm_beta_E);
end
%pensare se aggiungere sottoclasse o in generale un modo per salvare sta
%roba
if isnumeric(obj.beta)
    beta=@(x)obj.beta;
else
    beta=obj.beta;
end
forcing= @(x) obj.f(x) * dot(beta(x),el.gradPhi{j}(x)); %di nuovo cosa brutta

f=tau_E*quadrature(forcing);

end


