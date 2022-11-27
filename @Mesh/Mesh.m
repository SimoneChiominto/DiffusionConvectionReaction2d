classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here

    properties
        geom
        areaMax
        diamMax
        nDegreeOfFreedom
    end

    methods
        function obj = Mesh(Domain,BC,RefiningOptions)
            %MESH Construct an instance of this class
            %   Detailed explanation goes here
            addpath('@Mesh/Triangolatore/Long/bbtr30/')
            obj.geom = bbtr30(Domain,BC,RefiningOptions);

            % --------------------------------------------------
            % Rielaborazione dei prodotti del triangolatore
            % per un piu` agevole trattamento delle condizioni
            % al contorno
            % --------------------------------------------------

            obj.geom.elements.coordinates = obj.geom.elements.coordinates(...
                1:obj.geom.nelements.nVertexes,:);
            obj.geom.elements.triangles = obj.geom.elements.triangles(...
                1:obj.geom.nelements.nTriangles,:);
            obj.geom.elements.borders = obj.geom.elements.borders(...
                1:obj.geom.nelements.nBorders,:);
            obj.geom.elements.neighbourhood = obj.geom.elements.neighbourhood(...
                1:obj.geom.nelements.nTriangles,:);

            % --------------------------------------------------

            j  = 1;
            Dj = 1;
            for i=1:size(obj.geom.pivot.nodelist)
                if obj.geom.pivot.nodelist(i)==0
                    obj.geom.pivot.pivot(i)=j;
                    j = j+1;
                else
                    obj.geom.pivot.pivot(i)=-Dj;
                    Dj = Dj + 1;
                end
            end

            % --------------------------------------------------

            obj.geom.pivot.pivot = transpose(obj.geom.pivot.pivot);

            % --------------------------------------------------

            % geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
            % di Dirichlet e il corrispondente marker

            [X,I] = sort(obj.geom.pivot.Di(:,1));
            obj.geom.pivot.Di = obj.geom.pivot.Di(I,:);

            clear X I;

        end

        []=P2(obj) %declare method in external file

        function area = getAreaMax(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if isempty(obj.areaMax)
                area=max([obj.geom.support.TInfo.Area]);
                obj.areaMax=area;
                return
            end
            area=obj.areaMax;
            return
        end

        function h = getDiamMax(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if isempty(obj.diamMax)
                    x=obj.geom.elements.coordinates(obj.geom.elements.borders(:,1),1)-obj.geom.elements.coordinates(obj.geom.elements.borders(:,2),1);
                    y=obj.geom.elements.coordinates(obj.geom.elements.borders(:,1),2)-obj.geom.elements.coordinates(obj.geom.elements.borders(:,2),2);
                    h=max((x.^2+y.^2).^0.5);
                    obj.diamMax=h;
                return
            end
            h=obj.diamMax;
            return
        end

        function nDof = getDoF(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if isempty(obj.nDegreeOfFreedom)
                nDof=max(obj.geom.pivot.pivot);
                obj.nDegreeOfFreedom=nDof;                    
                return
            end
            nDof=obj.nDegreeOfFreedom;
            return
        end
    end
end

