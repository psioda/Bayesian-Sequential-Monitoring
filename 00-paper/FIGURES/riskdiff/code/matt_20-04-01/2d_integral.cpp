// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// P(a1 < X1 < b1, a2 < X2 < b2), (X1, X2) ~ N([0], [1   rho])
//                                            ([0], [rho   1])
class BiNormal: public MFunc
{
private:
    const double rho;
    double const1;  // 2 * (1 - rho^2)
    double const2;  // 1 / (2 * PI) / sqrt(1 - rho^2)
public:
    BiNormal(const double& rho_) : rho(rho_)
    {
        const1 = 2.0 * (1.0 - rho * rho);
        const2 = 1.0 / (2 * M_PI) / std::sqrt(1.0 - rho * rho);
    }

    // PDF of bivariate normal
    double operator()(Constvec& x)
    {
        double z = x[0] * x[0] - 2 * rho * x[0] * x[1] + x[1] * x[1];
        return const2 * std::exp(-z / const1);
    }
};

// [[Rcpp::export]]
Rcpp::List integrate_test2()
{
    BiNormal f(0.5);  // rho = 0.5
    Eigen::VectorXd lower(2);
    lower << -1, -1;
    Eigen::VectorXd upper(2);
    upper << 1, 1;
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return Rcpp::List::create(
        Rcpp::Named("approximate") = res,
        Rcpp::Named("error_estimate") = err_est,
        Rcpp::Named("error_code") = err_code
    );
}