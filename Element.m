classdef Element
    %ELEMENT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        coordinates
        B
        Binv
        phi
        gradPhi
        Fe
        FeInv
        Area
        nDoF
    end

    methods
        function obj = Element(coordinates,refElement)
            %ELEMENT Construct an instance of this class
            %   Detailed explanation goes here
            obj.coordinates = coordinates;
            dx(1) = coordinates(3,1)-coordinates(2,1);
            dx(2) = coordinates(1,1)-coordinates(3,1);
            dx(3) = coordinates(2,1)-coordinates(1,1);
            dy(1) = coordinates(2,2)-coordinates(3,2);
            dy(2) = coordinates(3,2)-coordinates(1,2);
            dy(3) = coordinates(1,2)-coordinates(2,2);

            obj.B=[dx(2), -dx(1); -dy(2) dy(1)]; 
            obj.Area=det(obj.B)/2;
            if obj.Area<=0
                error("coordinate non in senso antiorario")
            end
            obj.Binv =[dy(1),dx(1);dy(2),dx(2)]/(2*obj.Area);                                                            
            obj.Fe= @(x_hat) obj.B*x_hat+coordinates(3,:)'; 
            obj.FeInv = @(x) obj.Binv*(x-coordinates(3,:)');
            obj.nDoF=length(refElement.phi);
            for j=1:obj.nDoF
                obj.phi{j}=@(x)  refElement.phi{j}(obj.FeInv(x));
                obj.gradPhi{j}= @(x) obj.Binv'*refElement.gradPhi{j}(obj.FeInv(x));
            end

        end

    end
end

