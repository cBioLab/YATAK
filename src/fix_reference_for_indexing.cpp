#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
using namespace std;

char fixBase(char c) {
  const string alpha = "ACGT";
  for(int i = 0; i < 4; i++) {
    if(alpha[i] == c || tolower(alpha[i]) == c) {
      return alpha[i];
    }
  }
  return c;
}

int main(int argc, char **argv){

  if(argc < 2) {
    cerr << "Usage : " << endl;
    cerr << " ./fix_reference_for_indexing <reference.fasta> <fixed-reference.fasta>" << endl;
    return 0;
  }

  string referenceFilePath = string(argv[1]);
  string fixedReferenceFilePath = string(argv[2]);
  ifstream ifsRef(referenceFilePath);
  ofstream ofsFixRef(fixedReferenceFilePath);
  string referenceId, referenceString;
  string input;
  while(getline(ifsRef, input)) {
    if(input[0] == '>') {
      ofsFixRef << input << '\n';
      continue;
    }
    for(int i = 0; i < (int) input.size(); i++) {
      ofsFixRef << fixBase(input[i]);
    }
    ofsFixRef << '\n';
  }

  ofsFixRef.close();

  return 0;
}
