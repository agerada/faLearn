#include <Rcpp.h>
#include <map>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
List kmers(const CharacterVector& x, int kmer = 3) {
  // Simple kmer algorithm that returns a named R list with simple kmer strings
  // and counts
  std::string dna_string = as<std::string>(x); 
  std::map<std::string, unsigned long long int> kmer_dict; 
  unsigned long long int n = dna_string.size() + 1 - kmer;
  for(int i = 0; i < n; i++){
    auto kmer_i = dna_string.substr(i, kmer); 
    if(kmer_dict.find(kmer_i) == kmer_dict.end()){
      kmer_dict[kmer_i] = 1; 
    } else {
      kmer_dict[kmer_i] = kmer_dict[kmer_i] + 1; 
    }
  }
  std::vector<unsigned long long int> kmer_output;
  std::vector<std::string> kmer_string;
  for(auto i = kmer_dict.begin(), n = kmer_dict.end(); i != n; i++){
    kmer_string.push_back(i->first); 
    kmer_output.push_back(i->second); 
  }
  return List::create(Named("kmer_string") = kmer_string, 
                      Named("kmer_value") = kmer_output); 
}

// [[Rcpp::export]]
CharacterVector generate_kmer_perms(int k, const CharacterVector& bases){
  /* This is my own algorithm, that lexiconographically increments
   * a string of length k to find all permutations with repetition of 
   * provided bases. Returns an R CharacterVector of all permutations */
  std::string b = Rcpp::as<std::string>(bases); 
  std::sort(b.begin(), b.end()); 
  std::vector<std::string> output_s; 
  char max_char = b.back(); 
  char min_char = b.front(); 
  std::string min_string;  
  for(int i = 0; i < k; i++){ 
    min_string.push_back(min_char);
  }
  std::string max_string; 
  for(int i = 0; i < k; i++){
    max_string.push_back(max_char); 
  }
  
  output_s.push_back(min_string);
  
  while(output_s.back() < max_string){
    std::string working_string = output_s.back(); 
    auto working_char = working_string.rbegin(); 
    
    while(working_char != working_string.rend()){
      if(*working_char < max_char){
        break; 
      } else {
        *working_char = b.front(); 
        working_char++; 
      }
    }
    auto working_char_index = b.find(*working_char); 
    char new_char = b[working_char_index + 1]; 
    *working_char = new_char; 
    output_s.push_back(working_string); 
  }
  return Rcpp::wrap(output_s); 
  /*
  std::string s = Rcpp::as<std::string>(bases); 
  std::sort(s.begin(), s.end());
  std::vector<std::string> output_s; 
  do {
    output_s.push_back(s); 
  } while (std::next_permutation(s.begin(), s.end()));
  return Rcpp::wrap(output_s);
   */
}


std::map<std::string, unsigned long long int> generate_kmer_perm_dict(int k, std::string b = "ACTG") {
  // "Private" function to generate a map of all possible kmers initialised to 
  // counts of 0
  std::sort(b.begin(), b.end()); 
  std::map<std::string, unsigned long long int> output_s; 
  char max_char = b.back(); 
  char min_char = b.front(); 
  std::string min_string;  
  for(int i = 0; i < k; i++){ 
    min_string.push_back(min_char);
  }
  std::string max_string; 
  for(int i = 0; i < k; i++){
    max_string.push_back(max_char); 
  }
  
  output_s[min_string] = 0;

  while(output_s.rbegin()->first < max_string){
    std::string working_string = output_s.rbegin()->first; 
    auto working_char = working_string.rbegin(); 
    
    while(working_char != working_string.rend()){
      if(*working_char < max_char){
        break; 
      } else {
        *working_char = b.front(); 
        working_char++; 
      }
    }
    auto working_char_index = b.find(*working_char); 
    char new_char = b[working_char_index + 1]; 
    *working_char = new_char; 
    output_s[working_string] = 0; 
  }
  return output_s; 
}

List wrap_custom(const std::map<std::string, unsigned long long int> dict){ 
  std::vector<unsigned long long int> kmer_output;
  std::vector<std::string> kmer_string;
  for(auto i = dict.begin(), n = dict.end(); i != n; i++){
    kmer_string.push_back(i->first); 
    kmer_output.push_back(i->second); 
  }
  return List::create(Named("kmer_string") = kmer_string, 
                      Named("kmer_value") = kmer_output); 
}

List wrap_custom(const std::vector<unsigned long long int> v){
  return List::create(v);
}

std::map<std::string, unsigned long long int> make_kmer_paired_list(
    const std::string& x, int kmer, 
    std::map<std::string, unsigned long long int> kmer_dict = {}) {
  // "Private" function that generates a paired named R list of kmers and 
  // counts, functionally identical to kmers()
  //std::map<std::string, int> kmer_dict; 
  unsigned long long int n = x.size() + 1 - kmer;
  for(unsigned long long int i = 0; i < n; i++){
    auto kmer_i = x.substr(i, kmer); 
    if(kmer_dict.find(kmer_i) == kmer_dict.end()){
      kmer_dict[kmer_i] = 1; 
    } else if (kmer_dict[kmer_i] == 0){
      kmer_dict[kmer_i] = 1; 
    } else {
      kmer_dict[kmer_i] = kmer_dict[kmer_i] + 1; 
    }
  }
  return kmer_dict; 
}

// [[Rcpp::export]]
List kmers_pointed(const CharacterVector& x, int kmer = 3, 
                   bool simplify = false, 
                   bool anchor = true) {
  // "Public" function that returns an R list of kmers. By default this is anchored
  // with all possible kmers (if none recorded in genome then = 0). If anchor=false
  // then currently behaves identically to kmers()
  // Need to implement simplify = true, which should return just a vector of kmer
  // counts. 
  std::string dna_string = as<std::string>(x); 

  if (anchor) {
    std::map<std::string, unsigned long long int> mapped = generate_kmer_perm_dict(kmer, "ACTG"); 
    if (simplify) {
      std::map<std::string, unsigned long long int> temp_dict = make_kmer_paired_list(dna_string, kmer, mapped); 
      std::vector<unsigned long long int> kmer_output; 
      for(auto i = temp_dict.begin(), n = temp_dict.end(); i != n; i++){
        kmer_output.push_back(i->second);
      }
      return wrap_custom(kmer_output); 
    }
    else {
      return wrap_custom(make_kmer_paired_list(dna_string, kmer, mapped));
    }
  }
  else {
    if (simplify) {
      std::map<std::string, 
               unsigned long long int> temp_dict = make_kmer_paired_list(dna_string, kmer); 
      std::vector<unsigned long long int> kmer_output; 
      for (auto i = temp_dict.begin(), n = temp_dict.end(); i != n; i++){
        kmer_output.push_back(i->second); 
      }
      return wrap_custom(kmer_output); 
    }
    return wrap_custom(make_kmer_paired_list(dna_string, kmer)); 
  }
}
