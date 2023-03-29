#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
DataFrame Cpp_ModelMatrix(Formula formula, DataFrame df){
  Rcpp::Environment stats_env("package:stats");
  Rcpp::Function model_frame = stats_env["model.matrix.default"];
  DataFrame df_new = model_frame(Rcpp::_["object"] = formula, Rcpp::_["data"] = df);
  return df_new;
}
// DataFrame Cpp_ModelMatrix(Formula formula, DataFrame df){
//     Rcpp::Environment stats_env("package:stats");
//     Rcpp::Function model_frame = stats_env["model.frame"];
//     DataFrame df_new = model_frame(Rcpp::_["formula"] = formula, Rcpp::_["data"] = df);
//     return df_new;
// }


// https://stackoverflow.com/questions/22828361/rcpp-function-to-select-and-to-return-a-sub-dataframe
Function subset("[.data.frame");

// [[Rcpp::export]]
DataFrame SubsetDataFrame(DataFrame x, uvec y) { // If this errors change uvec -> IntegerVector and convert below...
  return subset(x, y, R_MissingArg);
}

// [[Rcpp::export]]
List SubsetAllData(const DataFrame x, const List idList){ //const vec& ids, const vec& assigns){
  vec ids = idList["id"], assign = idList["assign"];
  
  int n = ids.size(), n_a = assign.size();
  if(n != n_a){
    stop("'idList$id' is not of same length as 'idList$assign'.");
  }
  
  assign -= 1;                       // Indexes
  vec id_col = x["id"];              // The id column
  Rcpp::List out(n);                 // Store 
  for(int i = 0; i < n; i ++){
    int current = ids[i], current_assign = assign[i];
    uvec current_u = find(id_col == current) + 1; // Add one since uvec indexes from zero
    DataFrame ss = SubsetDataFrame(x, current_u);
    out[current_assign] = ss;
  }
  
  return out;
}

// [[Rcpp::export]]
DataFrame CreateMatrices(const Formula fk, Rcpp::DataFrame SubData, const std::string What, const std::string Resp){
  if(What != "Y"){
    return Cpp_ModelMatrix(fk, SubData);
  }else{
    vec y = SubData[Resp];
    return DataFrame::create(Named(Resp) = y);
  }
}

// List CreateFixedMatrices(const List& formulas, const DataFrame& df, const vec& ids){
//   int n = ids.size(), K = formulas.size();
//   // X will be a list of length n, each of which will be a sub-list of length K
//   vec df_id = df["id"];
//   // Initiate store.
//   List X = List(n);
//   for(int num = 0; num < n; num++){
//     // int current = *num;
//     int current = ids[num];
//     List Xi = List(K);
//     // Which ROWS in df belong to this id?
//     uvec this_id = find(df_id == ids.at(num));
//     DataFrame this_df = SubsetDataFrame(df, this_id);
//     for(int k = 0; k < K; k++){ // Subset for each k in 1,...,K.
//       Formula fk = formulas[k];
//       Xi[k] = Cpp_ModelMatrix(fk, this_df);
//     }
//     X[num] = Xi;
//   }
//   return X;
// }