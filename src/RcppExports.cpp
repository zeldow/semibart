// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// semibartcpp
List semibartcpp(arma::mat iX, arma::mat itrt, arma::vec iy, double sigma, int sigdf, double sigquant, double kfac, double power, double base, arma::vec meanb, double sigb, int ntree, int ndpost, arma::ivec inumcut, int iusequants, double binary_offset, int probitlink, int verbose, int printevery);
RcppExport SEXP semibart_semibartcpp(SEXP iXSEXP, SEXP itrtSEXP, SEXP iySEXP, SEXP sigmaSEXP, SEXP sigdfSEXP, SEXP sigquantSEXP, SEXP kfacSEXP, SEXP powerSEXP, SEXP baseSEXP, SEXP meanbSEXP, SEXP sigbSEXP, SEXP ntreeSEXP, SEXP ndpostSEXP, SEXP inumcutSEXP, SEXP iusequantsSEXP, SEXP binary_offsetSEXP, SEXP probitlinkSEXP, SEXP verboseSEXP, SEXP printeverySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type iX(iXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type itrt(itrtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type iy(iySEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type sigdf(sigdfSEXP);
    Rcpp::traits::input_parameter< double >::type sigquant(sigquantSEXP);
    Rcpp::traits::input_parameter< double >::type kfac(kfacSEXP);
    Rcpp::traits::input_parameter< double >::type power(powerSEXP);
    Rcpp::traits::input_parameter< double >::type base(baseSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type meanb(meanbSEXP);
    Rcpp::traits::input_parameter< double >::type sigb(sigbSEXP);
    Rcpp::traits::input_parameter< int >::type ntree(ntreeSEXP);
    Rcpp::traits::input_parameter< int >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type inumcut(inumcutSEXP);
    Rcpp::traits::input_parameter< int >::type iusequants(iusequantsSEXP);
    Rcpp::traits::input_parameter< double >::type binary_offset(binary_offsetSEXP);
    Rcpp::traits::input_parameter< int >::type probitlink(probitlinkSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type printevery(printeverySEXP);
    rcpp_result_gen = Rcpp::wrap(semibartcpp(iX, itrt, iy, sigma, sigdf, sigquant, kfac, power, base, meanb, sigb, ntree, ndpost, inumcut, iusequants, binary_offset, probitlink, verbose, printevery));
    return rcpp_result_gen;
END_RCPP
}
