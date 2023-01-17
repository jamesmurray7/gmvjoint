#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
DataFrame Cpp_ModelMatrix(Formula formula, DataFrame df){
    Rcpp::Environment stats_env("package:stats");
    Rcpp::Function model_frame = stats_env["model.frame"];
    DataFrame df_new = model_frame(Rcpp::_["formula"] = formula, Rcpp::_["data"] = df);
    return df_new;
}


// https://stackoverflow.com/questions/22828361/rcpp-function-to-select-and-to-return-a-sub-dataframe
Function subset("[.data.frame");

// [[Rcpp::export]]
DataFrame SubsetDataFrame(DataFrame x, uvec y) { // If this errors change uvec -> IntegerVector and convert below...
  return subset(x, y, R_MissingArg);
}

// [[Rcpp::export]]
List CreateFixedMatrices(const List& formulas, const DataFrame& df, const vec& ids){
  int n = ids.size(), K = formulas.size();
  // X will be a list of length n, each of which will be a sub-list of length K
  vec df_id = df["id"];
  // Initiate store.
  List X = List(n);
  for(int num = 0; num < n; num++){
    // int current = *num;
    int current = ids[num];
    List Xi = List(K);
    // Which ROWS in df belong to this id?
    uvec this_id = find(df_id == ids.at(num));
    DataFrame this_df = SubsetDataFrame(df, this_id);
    for(int k = 0; k < K; k++){ // Subset for each k in 1,...,K.
      Formula fk = formulas[k];
      Xi[k] = Cpp_ModelMatrix(fk, this_df);
    }
    X[num] = Xi;
  }
  return X;
}