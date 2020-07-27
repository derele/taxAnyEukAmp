#include <Rcpp.h>
using namespace Rcpp;

// From dada2 by Benjamin J Callahan
//------------------------------------------------------------------
//' Checks a vector of character sequences for whether they are entirely ACGT.
//'
//' @param seqs A \code{character} of candidate DNA sequences.
//' @return A \code{logical}. Whether or not each input character was ACGT only.
//' @export
// [[Rcpp::export]]
Rcpp::LogicalVector isACGT(std::vector<std::string> seqs) {
  unsigned int i, pos, strlen;
  bool justACGT;
  const char *cstr;
  Rcpp::LogicalVector isACGT(seqs.size());
  
  for(i=0; i<seqs.size(); i++) {
    justACGT = true;
    strlen = seqs[i].length();
    cstr = seqs[i].c_str();
    for(pos=0; pos<strlen; pos++) {
      if(!(cstr[pos] == 'A' || cstr[pos] == 'C' || cstr[pos] == 'G' || cstr[pos] == 'T')) {
        justACGT = false;
        break;
      }
    }
    isACGT(i) = justACGT;
  }
  return(isACGT);
}
