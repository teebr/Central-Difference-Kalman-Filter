classdef(Abstract) CDKF < matlab.System  & matlab.system.mixin.Propagates
    % Implement Central Difference Kalman Filter. This is designed only as a
    % superclass. To create a usable CDKF, create a new system object which
    % inherits this class.
    % In the subclass, overwrite the stateFcn and outputFcn methods with
    % the model equations for a specific application. dummy state and
    % output functions are defined in this class to use as a template for
    % to ensure the correct inputs and outputs.
    % The subclass can contain additional properties (e.g. values for
    % lookup tables used within the state/output functions). The setupImpl
    % and validatePropertiesImpl methods can be overwritten. However, the
    % first line of the new method should call the coreSetup/ coreValidate
    % method first, which initialises this class. 
    % The stepImpl, resetImpl and input/output specification methods should
    % not be changed.
       
    properties
        x0 = 0 %initial state vector
    end
    
    properties(Nontunable)
        Q = 1 %Process noise matrix
        R = 1 %Output noise matrix
        P0 = 1 %Initial covariance estimate matrix
        h = sqrt(3) %tuning parameter
        dt = 1 %time step       
    end

    
    properties (Nontunable, Access = protected)
        nStates %number of states, inferred from x0
        nOutputs %number of outputs, inferred from R
        nAug %total number of states including noise (2*nStates + 1)
        Wm %weighting matrix to calculate mean of sigma points
        Wc %weighting matrix to calculate covariance of sigma points
        xVec %indexing vector to get states from sigma points
        wVec %indexing vector to get process noise from sigma points
        vVec %indexing vector to get output noise from sigma points
        onesL %vector of ones, length equal to nAug
        onesL2%vector of ones, length equal to number of sigma points
    end
    
    properties (DiscreteState)
        %augmented sigma point matrix (each row is a state estimate, each
        %column is a set of sigma points
        Xa 
    end
       
    methods (Access = protected)
        %% MAIN METHODS:
        function validatePropertiesImpl(obj)
            coreValidate(obj)
        end
        
        function coreValidate(obj)
            if ~ismatrix(obj.Q) || (size(obj.Q,1) ~= size(obj.Q,2))
                error('Process noise Q must be a square matrix')
            end
            
            if ~ismatrix(obj.R) || (size(obj.R,1) ~= size(obj.R,2))
                error('Output noise R must be a square matrix')
            end
            
            %P can be left empty (becomes initialised based on Q)
            if isempty(obj.P0)
                warning('Initial Covariance P0 empty: will be initialised based on Q')
            elseif ~ismatrix(obj.P0) || (size(obj.P0,1) ~= size(obj.P0,2))
                error('Initial Covariance P0 must be a square matrix')
            end
            
            if ~iscolumn(obj.x0)
                error('Initial state estimate x0 must be a column vector')
            end
            
            if numel(obj.x0) ~= size(obj.Q,1)
                error('Number of states must match number of rows of Q')
            end
            
            if ~isscalar(obj.h) || obj.h < 0
                error('Tuning parameter h must be a scalar greater than 0')
            end
            
            if ~isscalar(obj.dt) || obj.dt < 0
                error('Time step dt must be a scalar greater than 0')
            end
            
        end
        
        function setupImpl(obj)
            %all the set-up for the main CDKF routine is in the parentSetup
            %method. This makes it easier to rewrite the setupImpl function
            %for the subclass in case new parameters have been added which
            %require setting up also.
            coreSetup(obj)
        end
        
        function coreSetup(obj)
            if isempty(obj.P0)
                obj.P0 = 10*obj.Q;
            end
            
            ns = numel(obj.x0); %number of states
            no = size(obj.R,1); %number of outputs

            % calculate coefficents and other useful things
            l  = 2*ns + no; %number of augmented states
            l2 = 2*l + 1; %number of sigma points
            wm = 1/(2*obj.h^2) * ones(l2,1);% weighting vector
            wm(1) = (obj.h^2 - l)/(obj.h^2);
            
            obj.nStates = ns;
            obj.nOutputs = no;
            obj.nAug = l;
            obj.Wm = wm; %mean weighting vector
            obj.Wc = diag(wm); % covariance weighing matrix
            obj.onesL = ones(1,l);
            obj.onesL2 = ones(1,l2);
            obj.xVec = 1:ns;
            obj.wVec = ns+1:l-no;
            obj.vVec = l-no+1:l;
        end
                    
        function [xp,yp] = stepImpl(obj,u,y)
            % when this function is called at step k, it performs two
            % distinct operations:
            %1)
            %   retrieve xn[k] from discrete states.
            %   y[k] is calculated using xn[k] and u[k]
            %   xp[k] updated using Kalman gain and outputted
            %2)
            %   xn[k+1] is calculated using xp[k] and u[k].
            %   xn[k+1] is stored as a discrete state, so that it can be
            %   retrieved at the next step where it is now xn[k]
            
            
            %_____________________MEASUREMENT UPDATE______________________
            %load in sigma points (calculated at the last step), and
            %calculate Kalman gain etc
            
            % state and output sigma points
            Xx = obj.Xa(obj.xVec,:);
            Xv = obj.Xa(obj.vVec,:);
            xn = Xx*obj.Wm; %mean of sigma points for time update
           
            Y = obj.outputFcn(Xx,u,Xv);
            yn = Y*obj.Wm;
            
            %calculate covariances and Kalman gain
            dx = Xx - xn*obj.onesL2;
            dy = Y - yn*obj.onesL2;
            Px = dx*obj.Wc*dx';
            Py = dy*obj.Wc*dy';
            Pxy = dx*obj.Wc*dy';
            
            Px = 0.5*(Px + Px'); %enforce symmetry
            Py = 0.5*(Py + Py'); %enforce symmetry
            K = Pxy / Py;
            
            %use Kalman gain to update states and covariances
            xp = xn + K*(y - yn); %corrected state estimate
            PxOld = Px;
            Px = Px - K*Pxy'; %or K*Py*K'
            Px = 0.5*(Px + Px'); %enforce symmetry
            yp = obj.outputFcn(xp,u,0);
            
            %______________________TIME UPDATE___________________________
            % Now we can calculate the sigma points for the next timestep.
            % Start by updating the augmented state and covariance
            % estimates (stored at discrete states)
            Pa = blkdiag(Px,obj.Q,obj.R); %total covariance

            %take square root of covariance
            [sa,chk] = chol(Pa,'lower');
            if chk > 0
%                 warning('Cholesky failed: using previous covariance')
                [sa,chk] = chol(blkdiag(PxOld,obj.Q,obj.R),'lower');
                if chk > 0
                    error('things are really not good here')
                end
            end
            
            %create sigma points, then decompose Xa into states and noises.
            hsa = obj.h*sa; %quicker to define on separate lines
            xa = [xp;zeros(obj.nStates+obj.nOutputs,1)];
            xao = xa*obj.onesL; %replicate xa in advance to improve execution speed

            XaOld = [xa,xao + hsa, xao - hsa]; %sigma points
            Xx = XaOld(obj.xVec,:); 
            Xw = XaOld(obj.wVec,:);
            Xv = XaOld(obj.vVec,:);
            XxNew = obj.stateFcn(Xx,u,Xw);
            obj.Xa = [XxNew;Xw;Xv];
        end
        
        %% Discrete State Methods
        function [sz,dt,cp] = getDiscreteStateSpecificationImpl(obj,name)
            n = numel(obj.x0); %this is called before setup function
            L = 2*n + 1;
            if strcmp(name,'Xa')
                sz = [L , 2*L + 1];
            else
                error('Unknown State')
            end
            dt = 'double';
            cp = false;
        end
        
        function resetImpl(obj)
            % Initialize discrete-state properties.
            
            ns = numel(obj.x0); %number of states
            no = size(obj.R,1); %number of outputs
            
            Pa = blkdiag(obj.P0,obj.Q,obj.R); %total covariance
            
            %take square root of covariance
            [sa,chk] = chol(Pa,'lower'); %Pa from
            if chk > 0
                error('Cholesky failed')
            end
            
            %create sigma points, then decompose Xa into states and noises.
            hsa = obj.h*sa; %quicker to do transpose and multiplication on separate lines
            
            ol = ones(1,2*numel(obj.x0) + 1); %number of sigma points
            xa = [obj.x0;zeros(ns+no,1)];
            xao = xa*ol; %replicate xa in advance to improve execution speed
            obj.Xa= [xa,xao + hsa, xao - hsa]; %sigma points
        end
        
        %% INPUT/ OUTPUT SPECIFICATION
        function validateInputsImpl(obj,u,y)
            % Validate inputs to the step method at initialization.
            if ~iscolumn(u)
                error('u must be a column vector')
            end
            
            if ~iscolumn(y)
                error('y must be a column vector')
            end
            
            if numel(y) ~= size(obj.R,1)
                error('Size of output noise R must match number of measurements')
            end
        end
        
        function num = getNumOutputsImpl(~)
            num = 2;
        end
        
        function [fz1,fz2] = isOutputFixedSizeImpl(~)
            fz1 = true;
            fz2 = true;
        end
        
        function [sz1 , sz2] = getOutputSizeImpl(obj)
            sz1 = size(obj.x0);
            sz2 = propagatedInputSize(obj,2);
        end
        
        function [dt1, dt2] = getOutputDataTypeImpl(~)
            dt1 = 'double';
            dt2 = 'double';
        end
        
        function [cp1,cp2] = isOutputComplexImpl(~)
            cp1 = false;
            cp2 = false;
        end

        function flag = supportsMultipleInstanceImpl(~)
            flag = true;
        end
        
        %% STATE UPDATE/ OUTPUT (NEED TO BE OVERWRITTEN BY SUBCLASS!)
        function x_new = stateFcn(obj,x,u,w)
            x_new = x*(1-obj.dt) + u + w;   
        end
        
        function y = outputFcn(obj,x,u,v)
            y = x-u*obj.nStates + v;
        end
        
    end
end
