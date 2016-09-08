classdef ManOptPoint < ManifoldPoint
% Returns a point on a generic manifold manifold
    properties
        manifold
    end
    
    methods
        % Constructor of the manifold point
        % In subclasses of this, an intial value should be provided and it
        % should be enforced to be on the manifold
        function manoptPt = ManOptPoint(mydims,manifold,x)
            manoptPt@ManifoldPoint(mydims)
            manoptPt.manifold = manifold;
            manoptPt.dims = mydims;
            assert(manifold.dim() == prod(mydims))
            if nargin < 3
                manoptPt.x = manoptPt.rand();
            else
                manoptPt.x = [];
            end
        end

        % Basic informational functions
        % Name to print for the manifold point
        function [nm] = name(self)
            nm = self.manifold.name();
        end

        % Operations on the point
        % The objects named u are tangent vectors represented simply by real
        % vectors
        % The objects named y are other ManifoldPoints

        % Inner product of d1 and d2 on this point's tangent space
        function [v] = inner(self,u1,u2)
            v = self.manifold(self.x,u1,u2);
        end

        % Norm of u on this point's tangent space
        function [v] = norm(self, u)
            v = self.manifold.norm(self.x,u);
        end

        % Manifold distance between this point and y
        function [v] = dist(self, y)
            v = self.manifold.norm(self.x,y.x);
        end

        % Typical distance between this point and other points on the manifold
        function [v] = typicaldist(self)
            v = self.manifold.typicaldist();
        end

        % Projects u (generic) to the tangent space of this point
        function [v] = proj(self, u)
            v = self.manifold.proj(self.x,u);
        end

        % Embeds the Euclidean gradient egrad to gradient on the manifold rgrad
        function [rgrad] = egrad2rgrad(self, egrad)
            rgrad = self.manifold.egrad2rgrad(self.x,egrad);
        end

        % Embeds the Euclidean hessian egrad to the hessian on the manifold rhess
        function [rhess] = ehess2rhess(self, ehess)
            rhess = self.manifold.ehess2rhess(self.x,ehess);
        end

        % Puts u in the tangent space of this point if off due to noise (similar to proj)
        function [v] = tangent(self, u)
            v = self.manifold.tanget(self.x,u);
        end
        
        % Takes the exponential based on this point along t*u in the tangent space
        function [v] = exp(self, u, t)
            v = self.manifold.exp(self.x,u,t);
        end

        % Cheaper version of the exponential
        function [v] = retr(self, u, t)
            v = self.manifold.retr(self.x,u,t);
        end

        % Takes log of this point towards y (tangent vector pointing towards y)
        function [v] = log(self, y)
            v = self.manifold.log(self.x,y.x);
        end

        % Overwrites this point with a random point on the manifold
        % TODO: Make a function that returns a 'gaussian like' distribution
        %       centered at this point
        function rand(self)
            self.x = self.manifold.rand();
        end

        % Returns a random unit length vector in the tangent space of this point
        function [v] = randvec(self)
            v = self.manifold.randvec(self.x);
        end

        % Returns linear combination of two points on the tangent space of this point
        function [v] = lincomb(self,a1,u1,a2,u2)
            v = self.manifold.lincomb(self.x,a1,u1,a2,u2);
        end

        % Returns the zero tangent vector of this point
        function [v] = zerovec(self)
            v = self.manifold.zerovec(self.x);
        end

        % Takes the transport of u in this point's tangent space into y's tangent space
        function [v] = transp(self,y,u)
            v = self.manifold.transp(self.x,y.x,u);
        end

        % Computes the geodesic mean of this point and y
        function [v] = pairmean(self, y)
            v = self.manifold.pairmean(self.x,y.x);
        end

        % Returns real column vector representation of this manifold point
        function [v] = vec(self)
            v = self.manifold.vec(self.x);
        end

        % Returns matrix representation of this point
        function [v] = mat(self)
            v = self.manifold.mat(self.x);
        end

    end

end
