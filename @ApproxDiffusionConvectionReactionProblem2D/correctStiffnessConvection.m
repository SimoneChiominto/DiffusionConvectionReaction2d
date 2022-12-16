function [diff,conv] = correctStiffnessConvection(obj,el,k,j)
%CONVECTIONCORRECTION Summary of this function goes here
%   Detailed explanation goes here
if isnumeric(obj.mu)
    mu=@(x) obj.mu;
else
    mu=obj.mu;
end

if isnumeric(obj.beta)
    beta=@(x) obj.beta;
else
    beta=obj.beta;
end


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

diffusion= @(x) (mu(x) .* trace(el.hessPhi{k}(x))) * dot(beta(x),el.gradPhi{j}(x));
diffusion_component=quadrature(diffusion);

convection = @(x) dot(beta(x), el.gradPhi{k}(x))* dot(beta(x),el.gradPhi{j}(x));%da controllare la formula
convection_component=quadrature(convection);

diff= tau_E*diffusion_component;
conv= tau_E*convection_component;

end

