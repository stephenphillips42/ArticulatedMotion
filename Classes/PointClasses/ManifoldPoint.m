classdef ManifoldPoint
% Returns a point on a generic manifold manifold
    properties
        dims % Dimension(s) of the manifold (a 1 by k vector)
        x % Actual representation of the point (depends on manifold)
    end
    
    methods
        % Constructor of the manifold point
        % In subclasses of this, an intial value should be provided and it
        % should be enforced to be on the manifold
        function manifoldPt = ManifoldPoint(mydims)
            manifoldPt.dims = mydims;
            manifoldPt.x = [];
        end

        % Basic informational functions
        % Name to print for the manifold point
        function [nm] = name(self)
            nm = sprintf('Manifold of dimension (%s)', num2str(self.dims));
        end
        
        % Print the total dimension of the manifold
        function [d] = dim(self)
            d = prod(self.dims);
        end

        % Returns a hash of this point for identification
        function [hsh] = hash(self)
            hsh = ['z' hashmd5(self.x(:))];
        end

        % Operations on the point
        % The objects named u are tangent vectors represented simply by real
        % vectors
        % The objects named y are other ManifoldPoints

        % Inner product of u1 and u2 on this point's tangent space
        function [v] = inner(self,d1,d2)
            error('Not Implemented')
        end

        % Norm of u on this point's tangent space
        function [v] = norm(self, u)
            error('Not implemented')
        end

        % Manifold distance between this point and y
        function [v] = dist(self, y)
            error('Not implemented')
        end

        % Typical distance between this point and other points on the manifold
        function [v] = typicaldist(self)
            error('Not implemented')
        end

        % Projects u (generic) to the tangent space of this point
        function [v] = proj(self, u)
            error('Not implemented')
        end

        % Embeds the Euclidean gradient egrad to gradient on the manifold rgrad
        function [rgrad] = egrad2rgrad(self, egrad)
            error('Not implemented')
        end

        % Embeds the Euclidean hessian egrad to the hessian on the manifold rhess
        function [rhess] = ehess2rhess(self, ehess)
            error('Not implemented')
        end

        % Puts u in the tangent space of this point if off due to noise (similar to proj)
        function [v] = tangent(self, u)
            error('Not implemented')
        end
        
        % Takes the exponential based on this point along t*u in the tangent space
        function [v] = exp(self, u, t)
            error('Not implemented')
        end

        % Cheaper version of the exponential
        function [v] = retr(self, u, t)
            error('Not implemented')
        end

        % Takes log of this point towards y (tangent vector pointing towards y)
        function [v] = log(self, y)
            error('Not implemented')
        end

        % Overwrites this point with a random point on the manifold
        % TODO: Make a function that returns a 'gaussian like' distribution
        %       centered at this point
        function rand(self)
            error('Not implemented')
        end

        % Returns a random unit length vector in the tangent space of this point
        function [v] = randvec(self)
            error('Not implemented')
        end

        % Returns linear combination of two points on the tangent space of this point
        function [v] = lincomb(self,a1,u1,a2,u2)
            error('Not implemented')
        end

        % Returns the zero tangent vector of this point
        function [v] = zerovec(self)
            error('Not implemented')
        end

        % Takes the transport of u in this point's tangent space into y's tangent space
        function [v] = transp(self,y,u)
            error('Not implemented')
        end

        % Computes the geodesic mean of this point and y
        function [v] = pairmean(self, y)
            error('Not implemented')
        end

        % Returns real column vector representation of this manifold point
        function [v] = vec(self)
            error('Not implemented')
        end

        % Returns matrix representation of this point
        function [v] = mat(self)
            error('Not implemented')
        end

    end

end
