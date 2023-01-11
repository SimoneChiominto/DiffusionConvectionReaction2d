function [] = generateLinearSystem(obj,options)
%GENERATELINEARSYSTEM Compute Stiffness matrix and forcing
%vector
arguments
    obj ApproxDiffusionConvectionReactionProblem2D
    options.massLumping         logical = true
    options.correctConvection   logical = true
end
%save useful variables
if isnumeric(obj.beta) && all(obj.beta==[0;0])
    options.correctConvection = false ;
end
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
%obj.stiffnessMatrix=spalloc(n_dof,n_dof,10*n_dof);
%Ad=spalloc(n_dof,n_d,10*n_dof);

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
    el=getEl(obj.mesh,obj.refElement,e);

    if any([isnumeric(mu),isnumeric(beta),isnumeric(gamma)]) && el.type=="P1"
        dx=getdx(el);
        dy=getdy(el);
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
                    
                if isnumeric(mu) && el.type~="P1"
                    diffusion= @(x) mu .* el.gradPhi{k}(x)'*el.gradPhi{j}(x);
                elseif isa(mu,"function_handle")
                    diffusion= @(x) mu(x) .* el.gradPhi{k}(x)'*el.gradPhi{j}(x);
                end

                if isnumeric(beta) && el.type~="P1"
                    convection = @(x) beta' * el.gradPhi{k}(x) * el.phi{j}(x);
                elseif isa(beta,"function_handle")
                    convection = @(x) beta(x)' * el.gradPhi{k}(x) * el.phi{j}(x);
                end

                if isnumeric(gamma) && el.type~="P1"
                    reaction= @(x) gamma * el.phi{k}(x) * el.phi{j}(x);
                elseif isa(gamma,"function_handle")
                    reaction= @(x) gamma(x) * el.phi{k}(x) * el.phi{j}(x);
                end

                %check if it is a degree of freedom
                if kk>0
                    if options.correctConvection
                        [diff_correction,conv_correction]=obj.correctStiffnessConvection(el,k,j);
                    else
                        diff_correction=0;
                        conv_correction=0;
                    end

                    if isnumeric(mu) && el.type=="P1"
                        diffusion_component= mu * (dy(k)*dy(j)+dx(k)*dx(j))/(4*el.Area);
                    else
                        diffusion_component=quadrature(diffusion);
                    end

                    if isnumeric(beta) && el.type=="P1"
                        convection_component= (beta(1)*dy(k)+beta(2)*dx(k))/6;
                    else
                        convection_component=quadrature(convection);
                    end

                    if isnumeric(gamma) && el.type=="P1"
                        reaction_component=(gamma*el.Area*(1+(j==k)))/12;
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

                    if isnumeric(mu) && el.type=="P1"
                        diffusion_component= mu * (dy(k)*dy(j)+dx(k)*dx(j))/(4*el.Area);
                    else
                        diffusion_component=quadrature(diffusion);
                    end

                    if isnumeric(beta) && el.type=="P1"
                        convection_component= (beta(1)*dy(k)+beta(2)*dx(k))/6;
                    else
                        convection_component=quadrature(convection);
                    end

                    if isnumeric(gamma) && el.type=="P1"
                        reaction_component=(gamma*el.Area*(1+(j==k)))/12;
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
            end %for k=1:3

            forcing_term= @(x) f(x)* el.phi{j}(x);
            if options.correctConvection
                f_correction=obj.correctForcingConvection(el,j);
            else
                f_correction=0;
            end
            b(jj) = b(jj)...
                +quadrature(forcing_term)+f_correction;


        end %if jj>0

    end %for j=1:3
end %for e= 1:obj.mesh.geom.nelements.nTriangles

%DA FARE (dovrebbe funzionare)
for e=1:length(obj.mesh.geom.pivot.Ne)
    %take the right boundary function
    g_n=boundaryFunctions{obj.mesh.geom.pivot.Ne(e,2)};

    l=obj.mesh.geom.pivot.Ne(e,1);
    i_b = obj.mesh.geom.elements.borders(l,1);
    i_e= obj.mesh.geom.elements.borders(l,2);
    ii_b=obj.mesh.geom.pivot.pivot(i_b);
    ii_e=obj.mesh.geom.pivot.pivot(i_e);

    x_b=obj.mesh.geom.elements.coordinates(i_b,:)';
    x_e=obj.mesh.geom.elements.coordinates(i_e,:)';
    
    if obj.mesh.elements{1}.type=="P2"
        i_m = obj.mesh.geom.elements.borders(l,5);
        ii_m=obj.mesh.geom.pivot.pivot(i_m);
        x_m=obj.mesh.geom.elements.coordinates(i_m,:)';

        if ii_b>0
            b_n(ii_b)=b_n(ii_b)+(2/15*g_n(x_b)+1/15*g_n(x_m)-1/30*g_n(x_e))...
                * norm([obj.mesh.geom.elements.coordinates(i_b,1)-obj.mesh.geom.elements.coordinates(i_e,1), obj.mesh.geom.elements.coordinates(i_b,2)-obj.mesh.geom.elements.coordinates(i_e,2)]);
        end
        if ii_e>0
            b_n(ii_e)=b_n(ii_e)+(2/15*g_n(x_e)+1/15*g_n(x_m)-1/30*g_n(x_b))...
                * norm([obj.mesh.geom.elements.coordinates(i_b,1)-obj.mesh.geom.elements.coordinates(i_e,1), obj.mesh.geom.elements.coordinates(i_b,2)-obj.mesh.geom.elements.coordinates(i_e,2)]);
        end
        b_n(ii_m)=b_n(ii_m)+(1/15*g_n(x_b)+8/15*g_n(x_m)+1/15*g_n(x_e))...
            * norm([obj.mesh.geom.elements.coordinates(i_b,1)-obj.mesh.geom.elements.coordinates(i_e,1), obj.mesh.geom.elements.coordinates(i_b,2)-obj.mesh.geom.elements.coordinates(i_e,2)]);


    else
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
    end
end %e=1:length(obj.mesh.geom.pivot.Ne)

u_d=zeros(length(obj.mesh.geom.pivot.Di),1);
for i=1:length(obj.mesh.geom.pivot.Di)
    %take the right boundary function
    g_d=boundaryFunctions{obj.mesh.geom.pivot.Di(i,2)};
    u_d(i)=g_d(obj.mesh.geom.elements.coordinates(obj.mesh.geom.pivot.Di(i,1),:)');
end

obj.stiffnessMatrix=sparse(JJ(1:count-1),KK(1:count-1),val(1:count-1),n_dof,n_dof);
Ad=sparse(JJ_D(1:countD-1),KK_D(1:countD-1),val_D(1:countD-1),n_dof,n_d);
%compute forcing vector
obj.forcingVector=b-Ad*u_d+b_n;

end %generateLinearSystem



% function obj = getEl(mesh,refElement,e)
% %if length(obj.elements)<e
%     coordinates=mesh.geom.elements.coordinates(mesh.geom.elements.triangles(e,1:3),:);
%     obj.coordinates = coordinates;
%     dx(1) = coordinates(3,1)-coordinates(2,1);
%     dx(2) = coordinates(1,1)-coordinates(3,1);
%     dx(3) = coordinates(2,1)-coordinates(1,1);
%     dy(1) = coordinates(2,2)-coordinates(3,2);
%     dy(2) = coordinates(3,2)-coordinates(1,2);
%     dy(3) = coordinates(1,2)-coordinates(2,2);
% 
%     obj.B=[dx(2), -dx(1); -dy(2) dy(1)];
%     obj.Area=det(obj.B)/2;
%     if obj.Area<=0
%         error("coordinate non in senso antiorario")
%     end
%     obj.Binv =[dy(1),dx(1);dy(2),dx(2)]/(2*obj.Area);
%     obj.Fe= @(x_hat) obj.B*x_hat+coordinates(3,:)';
%     obj.FeInv = @(x) obj.Binv*(x-coordinates(3,:)');
%     obj.nDoF=length(refElement.phi);
%     for j=1:obj.nDoF
%         obj.phi{j}=@(x)  refElement.phi{j}(obj.FeInv(x));
%         obj.gradPhi{j}= @(x) obj.Binv'*refElement.gradPhi{j}(obj.FeInv(x));
%         obj.hessPhi{j}= @(x) obj.B'*refElement.hessPhi{j}(x)*obj.B;
%     end
% 
%     switch obj.nDoF
%         case 3
%             obj.type="P1";
%         case 6
%             obj.type="P2";
%     end
%     %obj.elements(e)=el;
%     return
% %end
% %el=obj.elements(e);
% end
% 
%         function dx=getdx(obj)
%             dx(1) = obj.coordinates(3,1)-obj.coordinates(2,1);
%             dx(2) = obj.coordinates(1,1)-obj.coordinates(3,1);
%             dx(3) = obj.coordinates(2,1)-obj.coordinates(1,1);
%             return
%         end
% 
%         function dy=getdy(obj)
%             dy(1) = obj.coordinates(2,2)-obj.coordinates(3,2);
%             dy(2) = obj.coordinates(3,2)-obj.coordinates(1,2);
%             dy(3) = obj.coordinates(1,2)-obj.coordinates(2,2);
%             return
%         end
% 
