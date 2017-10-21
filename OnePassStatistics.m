classdef OnePassStatistics
    %This class can be used to compute Online Statistics (first 4 moments) for memory-consuming data
	%Data can be added sequentially and all current moments are available quickly through interface functions
    properties (Access = public)
        nElements
        M1
        M2
        M3
        M4
        MC
    end
    
    methods (Access = public)
        
        function obj = OnePassStatistics(dataDim)
            if nargin == 0
                dataDim = 1;
            end
            obj = obj.init(dataDim);
        end
        
        function obj = init(obj,dataDim)
            obj.nElements = 0;
            obj.M1 = zeros(dataDim);
            obj.M2 = zeros(dataDim);
            obj.M3 = zeros(dataDim);
            obj.M4 = zeros(dataDim);
            
            dataConsumpCovariance = prod(dataDim).^2 * 8 / 1024^3; %in GB
            
            if dataConsumpCovariance > 8 % check if bigger than 4 GB               
                warning(['Data dimension too big! No Covariance (approx. ' num2str(dataConsumpCovariance) 'GB) will be computed!']); %linearized storage for covariance
                obj.MC = NaN;
            else
                obj.MC  = zeros(prod(dataDim));
            end
            
        end
        
        function obj = addData(obj,val)
            if obj.nElements == 0
                obj = obj.init(size(val));
            end
            nNew = obj.nElements + 1;
            
            delta = val - obj.M1;
            delta_n = delta / nNew;
            %Update mean          
            obj.M1 = obj.M1 + delta_n;
            %update covariance
            if ~isnan(obj.MC)
                obj.MC = obj.MC + delta(:) * transpose(val(:) - obj.M1(:));
            end
            
            delta_nSq = delta_n.^2;
            term1 = obj.nElements * delta .* delta_n;
            
            %obj.M1 = obj.M1 +  delta_n;
            obj.M4 = obj.M4 +  term1 .* delta_nSq .* (nNew*nNew - 3*nNew + 3) + 6 * delta_nSq .* obj.M2 - 4 * delta_n .* obj.M3;
            obj.M3 = obj.M3 +  term1 .* delta_n .* (nNew - 2) - 3 .* delta_n .* obj.M2;
            obj.M2 = obj.M2 +  term1;
            
            
            obj.nElements = nNew;
        end
        
        %Overload numel function
        function out = numel(obj)
           out = obj.nElements;
        end
        
        %Mean
        function out = mean(obj)
            out = obj.M1;
        end
        
        %Variance
        function out = var(obj)
            out = obj.M2 ./ (obj.nElements);
        end
        
        %Covariance
        function out = cov(obj)
            out = obj.MC ./ (obj.nElements - 1);
        end
        
        %Standard deviation
        function out = std(obj)
            out = sqrt(var(obj));
        end
        
        %Skewness
        function out = skewness(obj)
            out = sqrt(obj.nElements) * obj.M3 ./ (obj.M2 .^ 1.5);
        end
        
        %Kurtosis
        function out = kurtosis(obj)
           out = obj.nElements * obj.M4 ./ (obj.M2.^2); 
           %out = out - 3; %Gaussian kurtosis normalization
        end
    end
    
end

