#include <Rcpp.h>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
using namespace Rcpp;

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

List wrap_custom(const std::map<unsigned long long int, unsigned long long int> dict){
  std::vector<unsigned long long int> kmer_index;
  std::vector<unsigned long long int> kmer_output;
  for(auto i = dict.begin(), n = dict.end(); i != n; i++){
    kmer_index.push_back(i->first);
    kmer_output.push_back(i->second);
  }
  return List::create(Named("kmer_index") = kmer_index,
                      Named("kmer_value") = kmer_output);
}

List wrap_custom(const std::vector<unsigned long long int> v){
  return List::create(v);
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

template<typename T>
bool is_valid_dna_string(T dna){
  return dna.size() > 0 ? true : false;
}

// [[Rcpp::export]]
std::string reverse_complement(std::string dna){

  std::string rev_comp;
  for(auto i = dna.rbegin(); i != dna.rend(); i++){
    switch(*i){
    case 'A':
      rev_comp.push_back('T');
      break;
    case 'T':
      rev_comp.push_back('A');
      break;
    case 'C':
      rev_comp.push_back('G');
      break;
    case 'G':
      rev_comp.push_back('C');
      break;
    default:
      rev_comp.push_back(*i);
      break;
    }
  }
  return rev_comp;
}

std::map<std::string, unsigned long long int> make_kmer_paired_list(
    const CharacterVector &x,
    int kmer,
    bool drop_n = false,
    bool canonical = true,
    std::map<std::string, unsigned long long int> kmer_dict = {}) {
  // "Private" function that generates a paired named R list of kmers and
  // counts, functionally identical to kmers()
  //std::map<std::string, int> kmer_dict;
  for (int i = 0; i < x.size(); i++) {
    std::string contig = as<std::string>(x[i]);

    if (!is_valid_dna_string(contig)) {
      throw std::range_error("Invalid DNA string (probably empty/NULL)");
    }
    unsigned long long int n = contig.size() + 1 - kmer;
    for(unsigned long long int j = 0; j < n; j++){
      auto kmer_i = contig.substr(j, kmer);

      // Check if kmer_i contains any non-ACTG characters
      if(drop_n){
        if(kmer_i.find_first_not_of("ACTG") != std::string::npos){
          continue;
        }
      }

      if(canonical){
        auto kmer_rev = reverse_complement(kmer_i);
        if(kmer_i > kmer_rev){
          kmer_i = kmer_rev;
        }
      }

      if(kmer_dict.find(kmer_i) == kmer_dict.end()){
        kmer_dict[kmer_i] = 1;
      } else if (kmer_dict[kmer_i] == 0){
        kmer_dict[kmer_i] = 1;
      } else {
        kmer_dict[kmer_i] = kmer_dict[kmer_i] + 1;
      }
    }
  }


  return kmer_dict;
}

std::map<unsigned long long int, unsigned long long int> convert_kmer_string_to_index(
    std::map<std::string, unsigned long long int> x,
    int k,
    int index) {
  std::map<unsigned long long int, unsigned long long int> output_dict;
  auto perms_dict = generate_kmer_perm_dict(k);
  for(auto i = perms_dict.begin(), n = perms_dict.end(); i != n; i++){
    i->second = index;
    index++;
  }
  for(auto i = x.begin(), n = x.end(); i != n; i++) {
    /* drop kmers that are not in permutations
     (i.e., ones that contain non-ACTG chars) */
    if(perms_dict.find(i->first) != perms_dict.end()) {
      unsigned long long int key_as_index = perms_dict[i->first];
      output_dict[key_as_index] = i->second;
    }
  }
  return output_dict;
}

//' Generates genome kmers
//'
//' @param x genome in string format
//' @param k kmer length
//' @param simplify returns a numeric vector of kmer counts, without associated string. This is useful to save memory, but should always be used with anchor = true.
//' @param canonical only record canonical kmers (i.e., the lexicographically smaller of a kmer and its reverse complement)
//' @param anchor includes unobserved kmers (with counts of 0). This is useful when generating a dense matrix where kmers of different genomes align.
//' @param clean_up only include valid bases (ACTG) in kmer counts (excludes non-coding results such as N)
//' @param key_as_int return kmer index (as "kmer_index") rather than the full kmer string. Useful for index-coded data structures such as libsvm.
//' @param starting_index the starting index, only used if key_as_int = TRUE.
//' @return list of kmer values, either as a list of a single vector (if simplify = TRUE), or as a named list containing "kmer_string" and "kmer_value".
//' @export
// [[Rcpp::export]]
List kmers(const CharacterVector& x,
          int k = 3,
          bool simplify = false,
          bool canonical = true,
          bool anchor = true,
          bool clean_up = true,
          bool key_as_int = false,
          bool starting_index = 1) {
 // "Public" function that returns an R list of kmers. By default this is anchored
 // with all possible kmers (if none recorded in genome then = 0). If anchor=false
 // then currently behaves identically to kmers()
 // If simplify = true, returns a numeric vector of kmer counts, without
 // associated string. This is useful to save memory, but should always be used
 // with anchor = true.
 // clean_up deals with missing data ("N") by dropping respective kmers
 // key_as_int converts the kmer string to an integer starting at starting_index, which
 // is useful for sparse matrices.
 // Note that if key_as_int=T, clean_up is implicit
 if (simplify & !anchor) warning("Simplifying but not anchoring - undefined behaviour");

 if (key_as_int) {
   auto string_key = make_kmer_paired_list(x, k, clean_up, canonical);
   auto int_key = convert_kmer_string_to_index(string_key, k, starting_index);
   return wrap_custom(int_key);
 }
 if (anchor) {
   std::map<std::string, unsigned long long int> mapped = generate_kmer_perm_dict(k, "ACTG");
   if (simplify) {
     std::map<std::string, unsigned long long int> temp_dict = make_kmer_paired_list(x, k, clean_up,
                                                                                     canonical,
                                                                                     mapped);
     std::vector<unsigned long long int> kmer_output;
     for(auto i = temp_dict.begin(), n = temp_dict.end(); i != n; i++){
       kmer_output.push_back(i->second);
     }
     return wrap_custom(kmer_output);
   }
   else {
     return wrap_custom(make_kmer_paired_list(x, k, clean_up, canonical, mapped));
   }
 }
 else {
   if (simplify) {
     std::map<std::string,
              unsigned long long int> temp_dict = make_kmer_paired_list(x, k, canonical, clean_up);
     std::vector<unsigned long long int> kmer_output;
     for (auto i = temp_dict.begin(), n = temp_dict.end(); i != n; i++){
       kmer_output.push_back(i->second);
     }
     return wrap_custom(kmer_output);
   }
   return wrap_custom(make_kmer_paired_list(x, k, clean_up));
 }
}

// [[Rcpp::export]]
bool kmers_to_libsvm(const CharacterVector &x,
                    const CharacterVector &target_path,
                    const CharacterVector &label = CharacterVector::create("0"),
                    int k = 3,
                    bool canonical = true) {
  //std::string dna_string = as<std::string>(x);
  std::ofstream file;
  std::string path = as<std::string>(target_path);
  file.open(path);
  file << as<std::string>(label) << " ";
  auto string_key = make_kmer_paired_list(x, k, true, canonical);
  auto int_key = convert_kmer_string_to_index(string_key, k, 1);
  for (auto const& i : int_key) {
    file << i.first << ":" << i.second << " ";
  }
  file << std::endl;
  file.close();
  return true;
}
