%> @file  spg_05262015.m
%> @brief Function that runs Spectral Gradient Projection Solver
%======================================================================
%> @brief It takes as function handles for the function and its gradient.
% It returns the optimal cost along with the optimal solution value
%> @param par initial value of variable
%> @param fn handle to function
%> @param gr handle to gradient
%> @param method optional scalar value to specify method
%> @param project optional projection function 
%> @param lower optional lower bound
%> @param upper optional upper bound
%> @param verbose integer if 1 outputs verbose messages
%> @retval sol structure containing all values at convergence
%> @b Author: 
%> Kayhan Batmanghelich
%>
%> @b Link: 
%> https://www.cbica.upenn.edu/sbia/software/
%> 
%> @b Contact: 
%> sbia-software@uphs.upenn.edu
%======================================================================
function  sol = spg_05262015(par, fn, gr, method, project, lower, upper,  verbose) 
  if (nargin<3)
	gr = @(x)x ;
	method = 3 ;
	project = @(x)x ;
	lower=-inf ;
	upper = inf ;
        warning('this part was not part of the original code') ;
  end
  %> control defaults
  ctrl.M=10;
  ctrl.maxit=1500;
  ctrl.gtol=1.e-05;
  ctrl.maxfeval=10000;
  ctrl.maximize= false;
  ctrl.trace = verbose;
  ctrl.triter=10;
  ctrl.eps=1e-7;
  ctrl.checkGrad_tol=1.e-06;
  ctrl.ftol = 1e-4 ;   % relative tolerance in f
  
  M	    = ctrl.M ;
  maxit    = ctrl.maxit ;
  gtol     = ctrl.gtol ;
  ftol     = ctrl.ftol ;
  maxfeval = ctrl.maxfeval ;
  maximize = ctrl.maximize ;
  trace    = ctrl.trace ;
  triter   = ctrl.triter ;
  eps      = ctrl.eps ;
  checkGrad_tol = ctrl.checkGrad_tol ; 

  if (any(isfinite(lower)) && length(lower)==1) 
	lower = repmat(lower, length(par),1) ;
  end
  if (any(isfinite(upper)) && length(upper)==1) 
	upper = repmat(upper, length(par),1) ;
  end
  


  grNULL = isempty(gr) ; 
  %fargs <- list(...)
  %#############################################
  if (~grNULL) 
    %grad_num <- diff(fn,par) ;  
    grad_analytic = feval(gr,par) ;
%    warning(' checking analytic gradient with numeric gradient is disable!') ; 
    %max_diff = max(abs((grad_analytic - grad_num) / (1 + abs(feval(par))))) ;
    %if(!max.diff < checkGrad.tol) {
    %  cat("Gradient check details:  max. relative difference in gradients= ",
    %	         max.diff,
    %	   "\n\n  analytic gradient:",  grad.analytic,
    %	   "\n\n  numerical gradient:", grad.num
    %	   )
    %  stop("Analytic gradient does not seem correct! See comparison above. ",
    %       "Fix it, remove it, or increase checkGrad.tol." )
    %  }
    %}
  end
  %> local function 
  %> Simple gr numerical approximation. Using func, f and eps from calling env.  	
  %> used when user does not specify gr.
  if (grNULL) gr  
	gr = @(x)defaultGradient(x,fn,eps) ;
  end


  %> This provides box constraints defined by upper and lower
  %> local functions defined only when user does not specify project.
  if (isempty(project)) 
       project = @(x)defaultProjection(x, lower, upper) ;
  end
  %#############################################

  %>  Initialization
  lmin = 1.e-10 ;
  lmax = 1.e10 ;
  iter = 0 ;  
  f_eval = 0 ;  
  geval  = 0 ;
  lastfv = repmat(-1.e99, M,1) ;
  fbest = NaN ;
 
  %> c() in next is for case of a 1x1 matrix value
  if (maximize)
  	func = -@(x)fn(x) ;
        grad = -@(x)gr(x) ;
  else
	func = @(x)fn(x) ;
   	grad = @(x)gr(x) ;
  end


  %> Project initial guess
  par2 = feval(project,par) ;

  if any(isnan(par2))
  	error('Failure in projecting intial guess!') ;
  else
      par = par2;
  end 
  
  pbest = par ;
 
  f = feval(func,par) ;

  f_eval = f_eval + 1 ;

  if ((~isnumeric(f)) || (1 ~= length(f)))
      error('function must return a scalar numeric value! : %d',f) ;
  elseif (isnan(f) || isinf(f) )
          error('Failure in initial function evaluation!') ;
  end
    
  f0 = f ; 
  fbest = f ; 
  f_rep = (-1)^maximize * f ;
 
  g = grad(par) ;
  
  geval = geval + 1 ;

  if any(isnan(g))
	error('Failure in intial gradient evaluation!') ;
  end 
 
  fbest =  f ;
  lastfv(1) = fbest  ;
 
  pg = feval(project,par - g) ;

  if any(isnan(pg))
	error('Failure in initial projection!') ;
  end 
 
  pg = pg - par ;

  pg2n = sqrt(sum(pg.*pg)) ;
  pginfn = max(abs(pg)) ;
  gbest = pg2n ;
  if (pginfn ~= 0) 
	lambda = min(lmax, max(lmin, 1/pginfn)) ;
  end
 
  if (trace) 
	fprintf(['iter: ',num2str(0), ' f-value: ', num2str(f0), ' pgrad: ',num2str(pginfn), '\n']) ;
  end

  %#######################
  %>  Main iterative loop
  %#######################
  lsflag = [] ;  %# for case when tol is already ok initially and while loop is skipped
  Obj_Hist = [f] ;
  fbest_Hist = [] ;
  delta_f = inf ;
  fbest_iter = 0 ;
  while( (pginfn > gtol) && (iter <= maxit) && (delta_f > ftol) ) 
      iter = iter + 1 ;
 
      d = feval(project,par - lambda * g) ;
 
      if any(isnan(d))   
        lsflag = 4 ;
        break
      end
 
      d = d - par ;
      gtd = sum(g .* d) ;
 
      if (isinf(gtd))
        lsflag = 4 ;
        break
      end
 
      nmls_ans = nmls(par, f, d, gtd, lastfv, f_eval , func, maxfeval) ;
      lsflag = nmls_ans.lsflag ;
 
      if (lsflag ~= 0) 
        break
      end
 
      f     = nmls_ans.f ;
      pnew  = nmls_ans.p ;
      f_eval = nmls_ans.f_eval ;
      lastfv(mod(iter,M) + 1) = f ;
      Obj_Hist = [Obj_Hist f] ;
 
      gnew = feval(grad,pnew) ;     
      geval = geval + 1 ;
 
      if any(isnan(gnew)) 
        lsflag = 3 ;
        break
      end
 
      s = pnew - par ;
      y = gnew - g ;
      sts = sum(s.*s) ;
      yty = sum(y.*y) ;
      sty = sum(s.*y) ;
 
      if (method==1) 
        if ((sts==0)  || (sty < 0))  
            lambda = lmax ;
        else 
            lambda = min(lmax, max(lmin, sts/sty)) ;
        end
      elseif (method==2) 
        if ((sty < 0) || (yty == 0)) 
            lambda = lmax ; 
        else 
            lambda = min(lmax, max(lmin, sty/yty)) ;
        end
      elseif (method==3) 
           if ((sts==0)  || (yty == 0)) 
                lambda = lmax ;
           else 
                lambda = min(lmax, max(lmin, sqrt(sts/yty))) ;
           end
      end
 
 
      par = pnew ;
      g   = gnew ;
 
      pg = feval(project,par - g) ;
 
      if any(isnan(pg))  
            lsflag =  4 ;
            break
      end

      pg = pg - par ;
      pg2n = sqrt(sum(pg.*pg)) ;
      pginfn = max(abs(pg)) ;
 
      f_rep = (-1)^maximize * f ;
      if (trace && (mod(iter,triter) == 0))
           fprintf(['iter: ',num2str(iter), ' f-value: ', num2str(f_rep), ' pgrad: ',num2str(pginfn), '\n']) ;
      end

      if (f < fbest) 
            fbest = f ;
            pbest = pnew ;
            gbest = pginfn ;
            fbest_Hist = [fbest_Hist fbest] ;
            fbest_iter = iter ;
            if (length(fbest_Hist)>10)
                delta_f = mean(abs(diff(fbest_Hist(end-10:end)))) ;
                delta_f = delta_f/fbest ;
                %if (delta_f < ftol)
                %    display('slow convergence !!!')
                %end
            end
      end
      
      %> last we update fbest was 20 iterations back! so let's check whether f is decreasing or not
      if ((iter - fbest_iter)>40) 
         delta_f = mean(abs(diff(Obj_Hist(end-20:end)))) ;
         delta_f = delta_f/f ; 
      end
      
  end      %> while condition loop concludes

  if (isempty(lsflag)) 
        warning('convergence tolerance satisified at intial parameter values.') ;
	    lsflag = 0 ;
  end
 
  conv1.type = nan ;
  conv1.message= 'Convergence status unknown' ;
  if (lsflag==0) 
    if (pginfn <= gtol) 
        conv1.type = 0 ;
        conv1.message= 'Successful convergence' ;
    end
    if (iter >= maxit)  
        conv1.type = 1 ;
        conv1. message = 'Maximum number of iterations exceeded' ;
    end
    if (delta_f <= ftol) 
        conv1.type = 6 ;
        conv1. message = 'Slow convergence ' ;
    end
  else 
      par = pbest ;
      f_rep = (-1)^maximize * fbest ;
      f = (-1)^maximize * fbest ;
      pginfn = gbest ;
      if (lsflag==1) 
        conv1.type = 3 ;
        conv1.message = 'Failure:  Error in function evaluation' ;
      end
      if (lsflag==2)
          conv1.type=2 ;
          conv1.message = 'Maximum function evals exceeded' ;
      end
      if (lsflag==3) 
            conv1.type=4 ;
            conv1.message = 'Failure:  Error in gradient evaluation' ;
      end
      if (lsflag==4) 
            conv1.type=5;
            conv1.message = 'Failure:  Error in projection' ;
      end
  end
 
  sol.par = par ; 
  sol.value = f_rep ; 
  sol.gradient = pginfn ; 
  sol.fn_reduction = (-1)^maximize * (f0 - f) ;
  sol.iter=iter ;
  sol.f_eval = f_eval ; 
  sol.convergence = conv1.type ;
  sol.message = conv1.message ;
  sol.Obj_Hist = Obj_Hist ;
end
 
 
 
% ################ local function
function   output = nmls(p, f, d, gtd, lastfv, f_eval, func, maxfeval)
    %> Non-monotone line search of Grippo with safe-guarded quadratic interpolation
    gamma = 1.e-04 ;
    fmax  = max(lastfv) ;
    alpha = 1 ;
    pnew = p + alpha*d ;
    fnew = feval(func, pnew) ;
    f_eval = f_eval + 1 ;
    %matlab_eps = 2.22e-16 ; 
    matlab_eps = realmin;
 
    if (isnan(fnew) )
        output.p = NaN ;
        output.f = NaN ;
        output.f_eval = NaN ;
        output.lsflag = 1 ;
    else
        output.lsflag = 0 ;
    end
 
    count=0;
    while(fnew > fmax + gamma*alpha*gtd && count<100) 
        if (alpha <= 0.1)
            alpha = alpha/2 ;
        else
            atemp = -(gtd*alpha^2) / (2*(fnew - f - alpha*gtd)) ;
            if ((atemp < 0.1) || (atemp > 0.9*alpha))
                atemp = alpha/2 ;
            end
            alpha = atemp ;
        end

    	pnew = p + alpha*d ;
        fnew =  feval(func, pnew) ;
    	f_eval = f_eval + 1 ;
        
        if (isnan(fnew) )
            output.p = NaN ;
            output.f = NaN ;
            output.f_eval = NaN ;
            output.lsflag = 1 ;
        end

        if (f_eval > maxfeval)
            output.p = NaN ;
            output.f = NaN ;
            output.f_eval = NaN ;
            output.lsflag = 2 ;
        end

        if (alpha < 10*matlab_eps)   % step is too small
            output.lsflag = 2 ;
            break ;
        end
            
    count=count+1;
    end  %while condition loop ends
%    fprintf(1,'%d\n',count);    
    output.p = pnew ;
    output.f = fnew ;
    output.f_eval = f_eval ;
    
end

%> this is default fucntion to compute the gradient in case it isnot provided
function df = defaultGradient(par,func,eps) 
    	df = zeros(length(par),1) ;
    	for i=1:length(par)
    	  dx = par ;
    	  dx(i) = dx(i) + eps ; 
    	  df(i) = (feval(func,dx) - feval(func,par))/eps ;
    	end
end

%> default projection which projects on the box constraints
function par = defaultProjection(par,lower,upper)
       %# Projecting to ensure that box-constraints are satisfied
       par(par < lower) = lower(par < lower) ;
       par(par > upper) = upper(par > upper) ;
end


