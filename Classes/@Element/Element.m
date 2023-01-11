classdef Element
    %ELEMENT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        coordinates
        B
        Binv
        phi
        gradPhi
        hessPhi
        Fe
        FeInv
        Area
        nDoF
        type
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
                obj.hessPhi{j}= @(x) obj.B'*refElement.hessPhi{j}(x)*obj.B;
            end

            switch obj.nDoF
                case 3
                    obj.type="P1";
                case 6
                    obj.type="P2";
            end

        end

        function h_E=getMaxLength(obj)
            len(1)=norm(obj.coordinates(3,:)-obj.coordinates(2,:));
            len(2)=norm(obj.coordinates(1,:)-obj.coordinates(3,:));
            len(3)=norm(obj.coordinates(2,:)-obj.coordinates(1,:));
            h_E=max(len);
            return
        end

        function dx=getdx(obj)
            dx(1) = obj.coordinates(3,1)-obj.coordinates(2,1);
            dx(2) = obj.coordinates(1,1)-obj.coordinates(3,1);
            dx(3) = obj.coordinates(2,1)-obj.coordinates(1,1);
            return
        end

        function dy=getdy(obj)
            dy(1) = obj.coordinates(2,2)-obj.coordinates(3,2);
            dy(2) = obj.coordinates(3,2)-obj.coordinates(1,2);
            dy(3) = obj.coordinates(1,2)-obj.coordinates(2,2);
            return
        end

    end
end

