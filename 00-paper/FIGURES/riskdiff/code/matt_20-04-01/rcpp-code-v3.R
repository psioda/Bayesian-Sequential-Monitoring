
require(Rcpp);
require(RcppArmadillo);
require(RcppNumerical);

sourceCpp(code='

	// [[Rcpp::depends(RcppArmadillo)]]
	// [[Rcpp::depends(RcppEigen)]]
	// [[Rcpp::depends(RcppNumerical)]]


	#include <RcppArmadillo.h>
	#include <RcppNumerical.h>   
	#include <vector>

	double absd (double z) { if (z<0) { return -1.0*z;} return z; }

	class dist_inner: public Numer::Func
	{
		private:
		
		double x;
	
    	double m0;
    	double s0;
		double h0;

		double h_y00,h_y01,h_y10,h_y11,h_lcmb0,h_lcmb1,a00,a01; 
		double c_y00,c_y01,c_y10,c_y11,c_lcmb0,c_lcmb1;
		
		double lognc;
	
		public:
    		dist_inner( double init_x, 	double init_m0, double init_s0, double init_h0, 
						double init_h_y00, double init_h_y01, double init_h_y10, double init_h_y11, double init_h_lcmb0, double init_h_lcmb1, double init_a00, double init_a01,
						double init_c_y00, double init_c_y01, double init_c_y10, double init_c_y11, double init_c_lcmb0, double init_c_lcmb1,
						double init_lognc) : 
						
                        x(init_x), m0(init_m0), s0(init_s0), h0(init_h0),
						h_y00(init_h_y00), h_y01(init_h_y01), h_y10(init_h_y10), h_y11(init_h_y11), h_lcmb0(init_h_lcmb0), h_lcmb1(init_h_lcmb1), a00(init_a00), a01(init_a01),
						c_y00(init_c_y00), c_y01(init_c_y01), c_y10(init_c_y10), c_y11(init_c_y11), c_lcmb0(init_c_lcmb0), c_lcmb1(init_c_lcmb1),
						lognc(init_lognc) {}

    	double operator()(const double & y) const
		{
			double f = -lognc;
			
			
			if (h0>0                             ) {f += -0.5*std::pow(absd((y-x)-m0)/s0,h0);              }
			if ((h_y00+h_y01)>0 and a00 > 0.00001) {f += a00*h_lcmb0 + a00*(h_y01*log(x) + h_y00*log(1-x));}
			if ((h_y10+h_y11)>0 and a01 > 0.00001) {f += a01*h_lcmb1 + a01*(h_y11*log(y) + h_y10*log(1-y));}			
			if ((c_y00+c_y01)>0                  ) {f +=    c_lcmb0 +     c_y01*log(x) + c_y00*log(1-x);   }
			if ((c_y10+c_y11)>0                  ) {f +=    c_lcmb1 +     c_y11*log(y) + c_y10*log(1-y);   }
			return exp( f );
    	}
	};   
   
	class dist_outer: public Numer::Func
	{
		private:
		
		double halfspace;		
		
    	double m0;
    	double s0;
		double h0;		
		
		double h_y00,h_y01,h_y10,h_y11,h_lcmb0,h_lcmb1,a00,a01; 
		double c_y00,c_y01,c_y10,c_y11,c_lcmb0,c_lcmb1;
		
		double lognc;
	


	
		public:
    		dist_outer( double init_halfspace, 	double init_m0, double init_s0, double init_h0, 
						double init_h_y00, double init_h_y01, double init_h_y10, double init_h_y11, double init_h_lcmb0, double init_h_lcmb1, double init_a00, double init_a01,
						double init_c_y00, double init_c_y01, double init_c_y10, double init_c_y11, double init_c_lcmb0, double init_c_lcmb1,
						double init_lognc) : 		
                        halfspace(init_halfspace), m0(init_m0), s0(init_s0), h0(init_h0),
						h_y00(init_h_y00), h_y01(init_h_y01), h_y10(init_h_y10), h_y11(init_h_y11), h_lcmb0(init_h_lcmb0), h_lcmb1(init_h_lcmb1), a00(init_a00), a01(init_a01),
						c_y00(init_c_y00), c_y01(init_c_y01), c_y10(init_c_y10), c_y11(init_c_y11), c_lcmb0(init_c_lcmb0), c_lcmb1(init_c_lcmb1),
						lognc(init_lognc) {}

		double operator()(const double& x) const
		{
			dist_inner f_inner(x,m0,s0,h0,  h_y00,h_y01,h_y10,h_y11,h_lcmb0,h_lcmb1,a00,a01,
											c_y00,c_y01,c_y10,c_y11,c_lcmb0,c_lcmb1,lognc);

			double lower=0.001,upper=0.999,err_est;
			int err_code;

			if (halfspace==1) { lower=x; }
			if (halfspace<0)  { lower=std::min(x - halfspace,upper); }
			if (lower==upper) { return 0;}
			return integrate(f_inner, lower, upper, err_est, err_code);
		}		
		
		
	};
   

	// [[Rcpp::export]]
	arma::vec integrate_bivar_dist(
			double m0, double s0, double h0, 
			double h_y00, double h_y01, double h_y10, double h_y11, double h_lcmb0, double h_lcmb1, double a00,double a01,
			double c_y00, double c_y01, double c_y10, double c_y11, double c_lcmb0, double c_lcmb1,double pp_calc=1)
	{
	
		double lower=0.001,upper=0.999;
		
		double err_est;
    	int err_code;	
	
		// compute normalizing constant;
		double h      = 0;
		double lognc  = 0;
		double pp     = 0;
		
		{
			dist_outer f_outer(h,m0,s0,h0,  
								h_y00,h_y01,h_y10,h_y11,h_lcmb0,h_lcmb1,a00,a01,
								c_y00,c_y01,c_y10,c_y11,c_lcmb0,c_lcmb1,lognc);
								
			lognc = log(integrate(f_outer, lower, upper, err_est, err_code));
			if (pp_calc==0)
			{
				arma::vec res(1);
				res[0] = lognc;
				return res;			
			}
		}
		
		if (pp_calc==1)
 		{
			h=1;
			dist_outer f_outer(h,m0,s0,h0,  
								h_y00,h_y01,h_y10,h_y11,h_lcmb0,h_lcmb1,a00,a01,
								c_y00,c_y01,c_y10,c_y11,c_lcmb0,c_lcmb1,lognc);
								
			pp = integrate(f_outer, lower, upper, err_est, err_code);											
		}
		else if (pp_calc<0)
 		{
			h=pp_calc;
			dist_outer f_outer(h,m0,s0,h0,  
								h_y00,h_y01,h_y10,h_y11,h_lcmb0,h_lcmb1,a00,a01,
								c_y00,c_y01,c_y10,c_y11,c_lcmb0,c_lcmb1,lognc);
								
			pp = integrate(f_outer, lower, upper, err_est, err_code);											
		}		
		
		
		arma::vec res(2);
		res[0] = pp;
		res[1] = lognc;
		return res;		
	}	

	double BetaBin(double y, double n, double alpha, double beta)
	{
		return exp(R::lchoose(n,y) + R::lbeta(y+alpha,n-y+beta) - R::lbeta(alpha,beta));
	}
	
	// [[Rcpp::export]]
	arma::vec compute_borrowing_parm(double rho0, double rho1,
			double  hy00, double  hy01, double  hy10, double  hy11, double h_lcmb0, double h_lcmb1,
			double   n0,  double   y01, double   n1, double   y11,  double c_lcmb0, double c_lcmb1)
	{	
		arma::vec borrow(3);
		
		double hn0 = hy00 + hy01;
		double hn1 = hy10 + hy11;		
		
		double alpha0   = hy01 + 0.5;
		double beta0    = hy00 + 0.5;
		
		double alpha1   = hy11 + 0.5;
		double beta1    = hy10 + 0.5;
		
		double obsProb  = BetaBin(y01,n0,alpha0,beta0)*BetaBin(y11,n1,alpha1,beta1);

		
		
		double c0  = 0;
		for(double x0=0;x0<=n0;x0++)
		{
			for(double x1=0;x1<=n1;x1++)
			{
				double p = BetaBin(x0,n0,alpha0,beta0)*BetaBin(x1,n1,alpha1,beta1);
				if (p<=obsProb) { c0+= p; }
			}
		}
		borrow[0]= c0;
		
		double maxborrow0 = std::min(1.0,rho0*n0/hn0);
		borrow[1] = std::min(1.0,c0*maxborrow0);

		double maxborrow1 = std::min(1.0,rho1*n1/hn1);
		borrow[2] = std::min(1.0,c0*maxborrow1);
	
		return borrow;
	}
	
	// [[Rcpp::export]]
	double compute_pos(double nMiss0, double nMiss1,double rho0, double rho1, double eCrit,
			double s_m0, double s_s0, double s_h0, 
			double  hy00, double  hy01, double  hy10, double  hy11, double h_lcmb0, double h_lcmb1, double a00,double a01,
			double   y00, double   y01, double   y10, double   y11, double c_lcmb0, double c_lcmb1)
	{
	
		double hn0 = hy00 + hy01;
		double hn1 = hy10 + hy11;
		
		double n0 = y00 + y01;
		double n1 = y10 + y11;		
	
		// compute predictive densities;
		
		arma::mat pDen(nMiss0+1,nMiss1+1);
		
		for (double yMiss01=0;yMiss01<=nMiss0;yMiss01++){ double yMiss00 = nMiss0-yMiss01;
		for (double yMiss11=0;yMiss11<=nMiss1;yMiss11++){ double yMiss10 = nMiss1-yMiss11;		
			
			double lcmb0 = R::lchoose(nMiss0,yMiss01) + c_lcmb0;
			double lcmb1 = R::lchoose(nMiss1,yMiss11) + c_lcmb1;
			
			double fy00 = y00 + yMiss00;
			double fy01 = y01 + yMiss01;
			double fy10 = y10 + yMiss10;
			double fy11 = y11 + yMiss11;			
			
			
			pDen(yMiss01,yMiss11) = integrate_bivar_dist(s_m0, s_s0, s_h0, 
						hy00, 		  hy01, 		 hy10, 			hy11,  h_lcmb0, h_lcmb1, a00, a01,
						fy00, 		  fy01, 		 fy10, 			fy11,    lcmb0,   lcmb1, 0)[0];
		}
		}
		
		pDen = pDen - pDen.max();
		pDen = exp(pDen);
		pDen = pDen / accu(pDen);
		
		
		
		double pos = 0;
		for (double yMiss01=0;yMiss01<=nMiss0;yMiss01++){ double yMiss00 = nMiss0-yMiss01;
		for (double yMiss11=0;yMiss11<=nMiss1;yMiss11++){ double yMiss10 = nMiss1-yMiss11;	

			double fy00 = y00 + yMiss00;
			double fy01 = y01 + yMiss01;
			double fy10 = y10 + yMiss10;
			double fy11 = y11 + yMiss11;
		
		
		
			// compute borrowing parameter;
				double alpha0   = hy01 + 0.5;
				double beta0    = hy00 + 0.5;
				double fn0      = n0+nMiss0;
				double obsProb0 = BetaBin(fy01,fn0,alpha0,beta0);
				double c00      = 0;
				for(double y=0;y<=fn0;y++)
				{
					double p = BetaBin(y,fn0,alpha0,beta0);
					if (p<=obsProb0) { c00+= p; }
				}
				
				double alpha1   = hy11 + 0.5;
				double beta1    = hy10 + 0.5;
				double fn1      = n1+nMiss1;
				double obsProb1 = BetaBin(fy11,fn1,alpha1,beta1);
				double c01      = 0;
				for(double y=0;y<=fn1;y++)
				{
					double p = BetaBin(y,fn1,alpha1,beta1);
					if (p<=obsProb1) { c01+= p; }
				}				

				double maxborrow0 = std::min(1.0,rho0*n0/hn0);
				double maxborrow1 = std::min(1.0,rho1*n1/hn1);
 				a00		        = std::min(1.0,c00*maxborrow0);
 				a01             = std::min(1.0,c01*maxborrow1);			
		
			double lcmb0 = R::lchoose(fn0,fy01);
			double lcmb1 = R::lchoose(fn1,fy11);		
		
			double pp = integrate_bivar_dist(s_m0, s_s0, s_h0, 
						hy00, hy01, hy10, hy11, h_lcmb0, h_lcmb1, a00, a01,
						fy00, fy01, fy10, fy11,   lcmb0,   lcmb1, 1)[0];
		
			if (pp>eCrit) { pos += pDen(yMiss01,yMiss11); }


		}
		}
		
		
		return pos;
	}	
	
')
