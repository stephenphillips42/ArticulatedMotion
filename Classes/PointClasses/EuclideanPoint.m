classdef EuclideanPoint < ManifoldPoint
% Returns a point on a generic manifold manifold. Most of this code from
% MANOPT
    
    methods
        % Constructor of the manifold point
        function euclidPt = EuclideanPoint(mydims,x)
            euclidPt@ManifoldPoint(mydims)
            if nargin < 2
                euclidPt.x = zeros(mydims);
            else
                assert(all(size(x) == mydims))
                euclidPt.x = x;
            end
        end

        % Typical distance between this point and other points on the manifold
        function [v] = typicaldist(self)
            v = sqrt(self.dim());
        end

        % Operations on the point
        % Inner product of d1 and d2 on this point's tangent space
        function [v] = inner(~,d1,d2)
            v = d1.'*d2;
        end

        % Norm of u on this point's tangent space
        function [v] = norm(~, u)
            v = norm(u);
        end

        % Manifold distance between this point and y
        function [v] = dist(self, y)
            v = norm(self.x - y.x);
        end


        % Projects u (generic) to the tangent space of this point
        function [v] = proj(~, u)
            v = u;
        end

        % Embeds the Euclidean gradient egrad to gradient on the manifold rgrad
        function [rgrad] = egrad2rgrad(~, egrad)
            rgrad = egrad;
        end

        % Embeds the Euclidean hessian egrad to the hessian on the manifold rhess
        function [rhess] = ehess2rhess(~, ehess)
            rhess = ehess;
        end

        % Puts u in the tangent space of this point if off due to noise (similar to proj)
        function [v] = tangent(~, u)
            v = u;
        end
        
        % Takes the exponential based on this point along t*u in the tangent space
        function [v] = exp(self, u, t)
            if nargin < 3
                v = self.x + u;
            else
                v = self.x + t*u;
            end
        end

        % Cheaper version of the exponential
        function [v] = retr(self, u, t)
            v = self.exp(u,t);
        end

        % Takes log of this point towards y (tangent vector pointing towards y)
        function [v] = log(self, y)
            v = self.x - y.x;
        end

        % Overwrites this point with a random point on the manifold
        function rand(self)
            self.x = randn(self.dims);
        end

        % Returns a random unit length vector in the tangent space of this point
        function [v] = randvec(self)
            v = normc(randn(self.dims));
        end

        % Returns linear combination of two points on the tangent space of this point
        function [v] = lincomb(~,a1,u1,a2,u2)
            v = a1*u1 + a2*u2;
        end

        % Returns the zero tangent vector of this point
        function [v] = zerovec(self)
            v = zeros(self.dims);
        end

        % Takes the transport of u in this point's tangent space into y's tangent space
        function [v] = transp(~,~,u)
            v = u;
        end

        % Computes the geodesic mean of this point and y
        function [v] = pairmean(self, y)
            v = EuclideanPoint(self.dims, 0.5*(self.x + y.x));
        end

        % Returns real column vector representation of this manifold point
        function [v] = vec(self)
            v = self.x(:);
        end

        % Returns matrix representation of this point
        function [v] = mat(self)
            v = self.x;
        end

    end

end
