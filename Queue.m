classdef Queue < handle
%--------------------------------------------------------------------------
% Class:        Queue < handle
%               
% Constructor:  Q = Queue(n);
%               Q = Queue(n,x0);
%               
% Properties:   (none)
%               
% Methods:                 Q.Enqueue(xi);
%               xi       = Q.Dequeue();
%               count    = Q.Count();
%               capacity = Q.Capacity();
%               bool     = Q.IsEmpty();
%               bool     = Q.IsFull();
%                          Q.Clear();
%               
% Description:  This class implements a first-in first-out (FIFO) queue of
%               numeric values
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 16, 2014
%--------------------------------------------------------------------------

    %
    % Private properties
    %
    properties (Access = private)
        k;                  % current number of elements
        n;                  % queue capacity
        x;                  % queue array
        tail;               % tail pointer
        head;               % head pointer
    end
    
    %
    % Public methods
    %
    methods (Access = public)
        %
        % Constructor
        %
        function this = Queue(n,x0)
            %----------------------- Constructor --------------------------
            % Syntax:       Q = Queue(n);
            %               Q = Queue(n,x0);
            %               
            % Inputs:       n is the maximum number of numbers that Q can
            %               hold
            %               
            %               x0 is a vector (of length <= n) of numbers to
            %               add to the queue during initialization
            %               
            % Description:  Creates a FIFO queue with capacity n
            %--------------------------------------------------------------
            
            % Initialize queue
            if (n == 0)
                Queue.ZeroCapacityError();
            end
            this.n = n;
            this.x = nan(n,1);
            
            if ((nargin == 2) && ~isempty(x0))
                % Insert given elements
                k0 = numel(x0);
                if (k0 > n)
                    % Queue overflow
                    Queue.OverflowError();
                else
                    this.x(1:k0) = x0(:);
                    this.SetLength(k0);
                    this.SetHead(1);
                    this.SetTail(k0 + 1);
                end
            else
                % Empty queue
                this.Clear();
            end
        end
        
        %
        % Add element to end of queue
        %
        function Enqueue(this,xi)
            %------------------------- Enqueue ----------------------------
            % Syntax:       Q.Enqueue(xi);
            %               
            % Inputs:       xi is a numeric value
            %               
            % Description:  Adds value xi to Q
            %--------------------------------------------------------------
            
            this.SetLength(this.k + 1);
            this.x(this.tail) = xi;
            this.SetTail(this.tail + 1);
        end
        
        %
        % Remove element from front of queue
        %
        function xi = Dequeue(this)
            %------------------------- Dequeue ----------------------------
            % Syntax:       xi = Q.Dequeue();
            %               
            % Outputs:      xi is a numeric value
            %               
            % Description:  Removes first value from Q
            %--------------------------------------------------------------
            
            this.SetLength(this.k - 1);
            xi = this.x(this.head);
            this.SetHead(this.head + 1);
        end
        
        %
        % Return number of elements in queue
        %
        function count = Count(this)
            %-------------------------- Count -----------------------------
            % Syntax:       count = Q.Count();
            %               
            % Outputs:      count is the number of values in Q
            %               
            % Description:  Returns the number of values in Q
            %--------------------------------------------------------------
            
            count = this.k;
        end
        
        %
        % Return queue capacity
        %
        function capacity = Capacity(this)
            %------------------------- Capacity ---------------------------
            % Syntax:       capacity = Q.Capacity();
            %               
            % Outputs:      capacity is the size of Q
            %               
            % Description:  Returns the maximum number of values that can 
            %               fit in Q
            %--------------------------------------------------------------
            
            capacity = this.n;
        end
        
        %
        % Check for empty queue
        %
        function bool = IsEmpty(this)
            %------------------------- IsEmpty ----------------------------
            % Syntax:       bool = Q.IsEmpty();
            %               
            % Outputs:      bool = {true,false}
            %               
            % Description:  Determines if Q is empty (i.e., contains zero
            %               values)
            %--------------------------------------------------------------
            
            if (this.k == 0)
                bool = true;
            else
                bool = false;
            end
        end
        
        %
        % Check for full queue
        %
        function bool = IsFull(this)
            %-------------------------- IsFull ----------------------------
            % Syntax:       bool = Q.IsFull();
            %               
            % Outputs:      bool = {true,false}
            %               
            % Description:  Determines if Q is full
            %--------------------------------------------------------------
            
            if (this.k == this.n)
                bool = true;
            else
                bool = false;
            end
        end
        
        %
        % Clear the queue
        %
        function Clear(this)
            %-------------------------- Clear -----------------------------
            % Syntax:       Q.Clear();
            %               
            % Description:  Removes all values from Q
            %--------------------------------------------------------------
            
            this.SetLength(0);
            this.SetTail(1);
            this.SetHead(1);
        end
    end
    
    %
    % Private methods
    %
    methods (Access = private)
        %
        % Set length
        %
        function SetLength(this,k)
            if (k < 0)
                Queue.UnderflowError();
            elseif (k > this.n)
                Queue.OverflowError();
            end
            this.k = k;
        end
        
        %
        % Set head pointer
        %
        function SetHead(this,head)
            % Wrap head pointer
            if (head > this.n)
                head = 1;
            end
            this.head = head;
        end
        
        %
        % Set tail pointer
        % 
        function SetTail(this,tail)
            % Wrap tail pointer
            if (tail > this.n)
                tail = 1;
            end
            this.tail = tail;
        end
    end
    
    %
    % Private static methods
    %
    methods (Access = private, Static = true)
        %
        % Overflow error
        %
        function OverflowError()
            error('Queue overflow');
        end
        
        %
        % Underflow error
        %
        function UnderflowError()
            error('Queue underflow');
        end
        
        %
        % No capacity error
        %
        function ZeroCapacityError()
            error('Queue with no capacity is not allowed');
        end
    end
end
