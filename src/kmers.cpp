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
//' Reverse complement of DNA string
//'
//' @param dna DNA string
//' @return reverse complement of DNA string
//' @export
//' @examples
//' reverse_complement("ATCG")
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
    unsigned int kmer,
    bool drop_n = false,
    bool canonical = true,
    bool squeeze = false,
    std::map<std::string, unsigned long long int> kmer_dict = {}) {
  // "Private" function that generates a paired named R list of kmers and
  // counts, functionally identical to kmers()
  //std::map<std::string, int> kmer_dict;
  for (int i = 0; i < x.size(); i++) {
    std::string contig = as<std::string>(x[i]);

    if (!is_valid_dna_string(contig)) {
      throw std::range_error("Invalid DNA string (probably empty/NULL)");
    }

    if (contig.size() < kmer) {
      continue;
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

  if (squeeze) {
    for (auto i = kmer_dict.begin(), n = kmer_dict.end(); i != n; i++){
      auto kmer_rev = reverse_complement(i->first);
      if (kmer_rev > i->first) {
        if (kmer_dict[kmer_rev] > 0) {
          throw std::range_error("A taxonomically larger reverse complement has recorded a count.");
        }
        kmer_dict.erase(kmer_rev);
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

std::map<unsigned long long int, unsigned long long int> convert_kmer_string_to_int_seq(
    std::map<std::string, unsigned long long int> x,
    unsigned int starting_index) {
  std::map<unsigned long long int, unsigned long long int> output_dict;
  for(auto i = x.begin(), n = x.end(); i != n; i++) {
    output_dict[starting_index] = i->second;
    starting_index++;
  }
  return output_dict;
}

//' Generates genome kmers
//'
//' @param x genome in string format
//' @param k kmer length
//' @param simplify returns a numeric vector of kmer counts,
//' without associated string. This is useful to save memory,
//' but should always be used with anchor = true.
//' @param canonical only record canonical kmers
//' (i.e., the lexicographically smaller of a kmer and its reverse complement)
//' @param squeeze remove non-canonical kmers
//' @param anchor includes unobserved kmers (with counts of 0).
//' This is useful when generating a dense matrix where kmers of different
//' genomes align.
//' @param clean_up only include valid bases (ACTG) in kmer counts
//' (excludes non-coding results such as N)
//' @param key_as_int return kmer index (as "kmer_index")
//' rather than the full kmer string. Useful for index-coded data structures
//' such as libsvm.
//' @param starting_index the starting index, only used if key_as_int = TRUE.
//' @return list of kmer values, either as a list of a single vector
//' (if simplify = TRUE), or as a named list containing "kmer_string" and
//' "kmer_value".
//' @export
//' @examples
//' kmers("ATCGCAGT")
// [[Rcpp::export]]
List kmers(const CharacterVector& x,
          int k = 3,
          bool simplify = false,
          bool canonical = true,
          bool squeeze = false,
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
                                                                                     squeeze,
                                                                                     mapped);
     std::vector<unsigned long long int> kmer_output;
     for(auto i = temp_dict.begin(), n = temp_dict.end(); i != n; i++){
       kmer_output.push_back(i->second);
     }
     return wrap_custom(kmer_output);
   }
   else {
     return wrap_custom(make_kmer_paired_list(x, k, clean_up, canonical, squeeze, mapped));
   }
 }
 else {
   if (simplify) {
     std::map<std::string,
              unsigned long long int> temp_dict = make_kmer_paired_list(x, k, clean_up, canonical, squeeze);
     std::vector<unsigned long long int> kmer_output;
     for (auto i = temp_dict.begin(), n = temp_dict.end(); i != n; i++){
       kmer_output.push_back(i->second);
     }
     return wrap_custom(kmer_output);
   }
   return wrap_custom(make_kmer_paired_list(x, k, clean_up));
 }
}

//' Converts a genome to kmers stored in libsvm format on disk
//'
//' @param x genome in string format
//' @param target_path path to store libsvm file (.txt)
//' @param label libsvm label
//' @param k kmer length
//' @param canonical only record canonical kmers
//' (i.e., the lexicographically smaller of a kmer and its reverse complement)
//' @param squeeze remove non-canonical kmers
//' @return boolean indicating success
//' @description
//' This function converts a single genome to a libsvm file containing kmer
//' counts. The libsvm format will be as follows:
//'
//' \preformatted{
//'   label 1:count 2:count 3:count ...
//' }
//' Label is optional and defaults to 0. The kmer counts are indexed by the
//' kmer index, which is the lexicographically sorted index of the kmer.
//' Libsvm is a sparse format.
//'
//' @seealso
//' For multiple genomes in a directory, processed in parallel, see [genomes_to_kmer_libsvm()]
//'
//' For more details on libsvm format, see
//' \url{https://xgboost.readthedocs.io/en/stable/tutorials/input_format.html}
//'
//' @export
//' @examples
//' temp_libsvm_path <- tempfile(fileext = ".txt")
//' genome_to_libsvm("ATCGCAGT", temp_libsvm_path)
//' readLines(temp_libsvm_path)
// [[Rcpp::export]]
bool genome_to_libsvm(const CharacterVector &x,
                    const CharacterVector &target_path,
                    const CharacterVector &label = CharacterVector::create("0"),
                    int k = 3,
                    bool canonical = true,
                    bool squeeze = false) {
  std::ofstream file;
  std::string path = as<std::string>(target_path);
  file.open(path);
  file << as<std::string>(label) << " ";
  std::map<unsigned long long int, unsigned long long int> int_key;
  if (!squeeze) {
    auto string_key = make_kmer_paired_list(x, k, true, canonical, squeeze);
    auto int_key = convert_kmer_string_to_index(string_key, k, 1);
    for (auto const& i : int_key) {
      file << i.first << ":" << i.second << " ";
    }
  } else {
    // If we're squeezing, we need a bit of extra work.
    // Provide a kmer_perms dict as input to ensure that we have the correct
    // alignment. Then just convert the str component (map->first) to a sequence
    // of integers. We do not use the convert_kmer_string_to_index function
    // since this maps to the full kmer set, which is not what we want.
    auto perms_dict = generate_kmer_perm_dict(k);
    auto string_key = make_kmer_paired_list(x, k, true, canonical, squeeze, perms_dict);
    auto int_seq_key = convert_kmer_string_to_int_seq(string_key, 1);

    // We now also need to remove 0-count kmers
    for (auto it = int_seq_key.cbegin(); it != int_seq_key.cend();) {
      if (it->second == 0) {
        it = int_seq_key.erase(it++);
      } else {
        ++it;
      }
    }
    for (auto const& i : int_seq_key) {
      file << i.first << ":" << i.second << " ";
    }
  }
  file << std::endl;
  file.close();
  return true;
}

std::vector<std::string> make_squeezed_mers(int k) {
  auto perms_dict = generate_kmer_perm_dict(k);
  for (auto i = perms_dict.begin(), n = perms_dict.end(); i != n; i++){
    auto kmer_rev = reverse_complement(i->first);
    if (kmer_rev > i->first) {
      perms_dict.erase(kmer_rev);
    }
  }
  std::vector<std::string> output;
  for (auto i = perms_dict.begin(), n = perms_dict.end(); i != n; i++){
    output.push_back(i->first);
  }
  return output;
}

std::vector<std::string> make_unsqueezed_mers(int k) {
  auto perms_dict = generate_kmer_perm_dict(k);
  std::vector<std::string> output;
  for (auto i = perms_dict.begin(), n = perms_dict.end(); i != n; i++){
    output.push_back(i->first);
  }
  return output;
}

//' Generates all permutations of squeezed kmers
//' @param k kmer length
//' @return vector of squeezed kmers
//' @export
//' @examples
//' squeezed_mers(3)
// [[Rcpp::export]]
StringVector squeezed_mers(int k = 3) {
  auto output = make_squeezed_mers(k);
  StringVector output_sv(output.size());
  output_sv = output;
  return output_sv;
}

//' Generates all permutations of unsqueezed kmers
//' @param k kmer length
//' @return vector of unsqueezed kmers
//' @export
//' @examples
//' unsqueezed_mers(3)
// [[Rcpp::export]]
StringVector unsqueezed_mers(int k = 3) {
  auto output = make_unsqueezed_mers(k);
  StringVector output_sv(output.size());
  output_sv = output;
  return output_sv;
}

//' Get str conversion of squeezed kmer using index
//' @param x integer vector of kmer indices
//' @param k kmer length
//' @param starting_index starting index (libsvm is usually indexed starting at 1)
//' @return vector of squeezed kmer strings
//' @export
//' @examples
//' squeezed_index_to_str(2, k = 3)
// [[Rcpp::export]]
StringVector squeezed_index_to_str(IntegerVector x,
                               int k,
                               unsigned int starting_index = 1) {
  auto squeezed_mers = make_squeezed_mers(k);
  std::vector<std::string> output;
  for (auto i = x.begin(); i != x.end(); i++) {
    if ((*i) - starting_index >= squeezed_mers.size()) {
      throw std::range_error("Index out of bounds");
    }
    output.push_back(squeezed_mers[(*i) - starting_index]);
  }
  StringVector output_sv(output.size());
  output_sv = output;
  return output_sv;
}

//' Get str conversion of unsqueezed kmer using index
//' @param x integer vector of kmer indices
//' @param k kmer length
//' @param starting_index starting index (libsvm is usually indexed starting at 1)
//' @return vector of unsqueezed kmer strings
//' @export
//' @examples
//' unsqueezed_index_to_str(2, k = 3)
// [[Rcpp::export]]
StringVector unsqueezed_index_to_str(IntegerVector x,
                               int k,
                               unsigned int starting_index = 1) {
  auto unsqueezed_mers = make_unsqueezed_mers(k);
  std::vector<std::string> output;
  for (auto i = x.begin(); i != x.end(); i++) {
    if ((*i) - starting_index >= unsqueezed_mers.size()) {
      throw std::range_error("Index out of bounds");
    }
    output.push_back(unsqueezed_mers[(*i) - starting_index]);
  }
  StringVector output_sv(output.size());
  output_sv = output;
  return output_sv;
}
