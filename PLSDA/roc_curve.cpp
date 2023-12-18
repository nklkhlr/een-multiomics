#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;


float min(NumericVector x) {
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  return *it;
}


float max(NumericVector x) {
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  return *it;
}


// [[Rcpp::export]]
NumericMatrix roc_curve(
	NumericVector x, IntegerVector y,
	int n_points
) {
	int tp, fp, tn, fn;
	NumericMatrix roc(n_points + 1, 2);
	float current_threshold = max(x) + .01;
	// ensuring that we always go from 0 => 1 regardless of prediction values
	float step_size = (current_threshold - (min(x) + .01))/n_points;
	int i, j;
	for (i = 0; i < n_points; i++) {
		tp = 0;
		fp = 0;
		tn = 0;
		fn = 0;
		for (j = 0; j < x.size(); j++) {
			if (x[j] >= current_threshold) {
				if (y[j] == 1) tp += 1;
				else fp += 1;
			} else {
				if (y[j] == 2) 	tn += 1;
				else fn += 1;
			}
		}
		// true positive rate
		roc(i, 1) = float(tp)/(tp + fn);
		// false positive rate
		roc(i, 0) = float(fp)/(fp + tn);
		current_threshold -= step_size;
	}
	roc(n_points, 1) = 1;
	roc(n_points, 0) = 1;
	
	return roc;
}
