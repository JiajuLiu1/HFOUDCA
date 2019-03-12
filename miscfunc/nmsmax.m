function [x, fmax, nf] = nmsmax(fun, x, stopit, savit, varargin)
x0 = x(:); 
n = length(x0);
if nargin < 3 | isempty(stopit), stopit(1) = 1e-3; end
tol = stopit(1);  % Tolerance for cgce test based on relative size of simplex.
if length(stopit) == 1, stopit(2) = inf; end  % Max no. of f-evaluations.
if length(stopit) == 2, stopit(3) = inf; end  % Default target for f-values.
if length(stopit) == 3, stopit(4) = 0; end    % Default initial simplex.
if length(stopit) == 4, stopit(5) = 0; end    % Default: show progress.
trace  = stopit(5);
if nargin < 4, savit = []; end                   % File name for snapshots.

V = [zeros(n,1) eye(n)];
f = zeros(n+1,1);
V(:,1) = x0; f(1) = feval(fun,x,varargin{:});
fmax_old = f(1);

if trace, fprintf('f(x0) = %9.4e\n', f(1)), end

k = 0; m = 0;

% Set up initial simplex.
scale = max(norm(x0,inf),1);
if stopit(4) == 0
   % Regular simplex - all edges have same length.
   % Generated from construction given in reference [18, pp. 80-81] of [1].
   alpha = scale / (n*sqrt(2)) * [ sqrt(n+1)-1+n  sqrt(n+1)-1 ];
   V(:,2:n+1) = (x0 + alpha(2)*ones(n,1)) * ones(1,n);
   for j=2:n+1
       V(j-1,j) = x0(j-1) + alpha(1);
       x(:) = V(:,j); f(j) = feval(fun,x,varargin{:});
   end
else
   % Right-angled simplex based on co-ordinate axes.
   alpha = scale*ones(n+1,1);
   for j=2:n+1
       V(:,j) = x0 + alpha(j)*V(:,j);
       x(:) = V(:,j); f(j) = feval(fun,x,varargin{:});
   end
end
nf = n+1;
how = 'initial  ';

[temp,j] = sort(f);
j = j(n+1:-1:1);
f = f(j); V = V(:,j);

alpha = 1;  beta = 1/2;  gamma = 2;

while 1    
k = k+1;

    fmax = f(1);
    if fmax > fmax_old
       if ~isempty(savit)
          x(:) = V(:,1); eval(['save ' savit ' x fmax nf'])
       end
       if trace
          fprintf('Iter. %2.0f,', k)
          fprintf(['  how = ' how '  ']);
          fprintf('nf = %3.0f,  f = %9.4e  (%2.1f%%)\n', nf, fmax, ...
                  100*(fmax-fmax_old)/(abs(fmax_old)+eps))
       end
    end
    fmax_old = fmax;

    %%% Three stopping tests from MDSMAX.M

    % Stopping Test 1 - f reached target value?
    if fmax >= stopit(3)
       msg = ['Exceeded target...quitting\n'];
       break  % Quit.
    end

    % Stopping Test 2 - too many f-evals?
    if nf >= stopit(2)
       msg = ['Max no. of function evaluations exceeded...quitting\n'];
       break  % Quit.
    end

    % Stopping Test 3 - converged?   This is test (4.3) in [1].
    v1 = V(:,1);
    size_simplex = norm(V(:,2:n+1)-v1(:,ones(1,n)),1) / max(1, norm(v1,1));
    if size_simplex <= tol
       msg = sprintf('Simplex size %9.4e <= %9.4e...quitting\n', ...
                      size_simplex, tol);
       break  % Quit.
    end

    %  One step of the Nelder-Mead simplex algorithm
    %  NJH: Altered function calls and changed CNT to NF.
    %       Changed each `fr < f(1)' type test to `>' for maximization
    %       and re-ordered function values after sort.

    vbar = (sum(V(:,1:n)')/n)';  % Mean value
    vr = (1 + alpha)*vbar - alpha*V(:,n+1); x(:) = vr; fr = feval(fun,x,varargin{:});
    nf = nf + 1;
    vk = vr;  fk = fr; how = 'reflect, ';
    if fr > f(n)
            if fr > f(1)
               ve = gamma*vr + (1-gamma)*vbar; x(:) = ve; fe = feval(fun,x,varargin{:});
               nf = nf + 1;
               if fe > f(1)
                  vk = ve; fk = fe;
                  how = 'expand,  ';
               end
            end
    else
            vt = V(:,n+1); ft = f(n+1);
            if fr > ft
               vt = vr;  ft = fr;
            end
            vc = beta*vt + (1-beta)*vbar; x(:) = vc; fc = feval(fun,x,varargin{:});
            nf = nf + 1;
            if fc > f(n)
               vk = vc; fk = fc;
               how = 'contract,';
            else
               for j = 2:n
                   V(:,j) = (V(:,1) + V(:,j))/2;
                   x(:) = V(:,j); f(j) = feval(fun,x,varargin{:});
               end
               nf = nf + n-1;
               vk = (V(:,1) + V(:,n+1))/2; x(:) = vk; fk = feval(fun,x,varargin{:});
               nf = nf + 1;
               how = 'shrink,  ';
            end
    end
    V(:,n+1) = vk;
    f(n+1) = fk;
    [temp,j] = sort(f);
    j = j(n+1:-1:1);
    f = f(j); V = V(:,j);

end   %%%%%% End of outer (and only) loop.

% Finished.
if trace, fprintf(msg), end
x(:) = V(:,1);