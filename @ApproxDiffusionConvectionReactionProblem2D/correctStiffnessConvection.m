function [diff,conv] = correctStiffnessConvection(obj,el,k,j)
%CONVECTIONCORRECTION Summary of this function goes here
%   Detailed explanation goes here


if length(obj.mesh.geom.elements.triangles(1,:))==3
    m_k = 1/3;
elseif length(obj.mesh.geom.elements.triangles(1,:))==6
    m_k = 1/24;
end

quadrature= @(f) quadrature_ref(@(x_hat) f(el.Fe(x_hat)) )*2*el.Area;

h_E=el.getMaxLength();
norm_beta_E = quadrature(@(x) norm(obj.beta(x))^2) .^0.5 ;
norm_mu_E = quadrature(@(x) obj.mu(x).^2).^0.5;
Pe_h = m_k*norm_beta_E*h_E / (2*norm_mu_E);
if Pe_h<=1
    tau_E=m_k*(h_E^2)/(4*norm_mu_E);
else
    tau_E=h_E/(2*norm_beta_E);
end
%pensare se aggiungere sottoclasse o in generale un modo per salvare sta
%roba
diffusion= @(x) (obj.mu(x) .* trace(el.hessPhi{k}(x))) * dot(obj.beta(x),el.gradPhi{j}(x));
convection = @(x) dot(obj.beta(x), el.gradPhi{k}(x))* dot(obj.beta(x),el.gradPhi{j}(x));%da controllare la formula

diff= tau_E*quadrature(diffusion);
conv= tau_E*quadrature(convection);


end

