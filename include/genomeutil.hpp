#ifndef INCLUDED_GENOMEUTIL
#define INCLUDED_GENOMEUTIL

#include <string>

int getNumberFromBase(char base) {
  const std::string alphabet = "ACGT";
  for(int i = 0; i < 4; i++) {
    if(alphabet[i] == base) return i;
  }
  return -1;
}

char getRCBase(char base) {
  switch(base) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default : return 'N';
  }
}

bool isValidRead(const std::string& read) {
  const std::string alphabet = "ACGT";
  int readlen = read.size();
  for(int i = 0; i < readlen; i++){
    if(std::find(alphabet.begin(), alphabet.end(), read[i]) == alphabet.end()) return false;
  }
  return true;
}

std::string getRCRead(const std::string& read) {
  int readlen = read.size();
  std::string result;
  for(int i = 0; i < readlen; i++) {
    result += getRCBase(read[i]);
  }
  reverse(result.begin(), result.end());
  return result;
}

#endif