#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sdsl/cst_sct3.hpp>

#include <genomeutil.hpp>

//// namespace
using namespace sdsl;
using namespace std;
using namespace std::chrono;

//// type definition
using timer = std::chrono::high_resolution_clock;
using CST = cst_sct3<>;

using pos_type = int64_t;

//// constant
const string alphabet = "ACGT";
const char termCharacter = '$';

void construct_refbv(CST& cst, pos_type refSize, const string& refbvFilePath) {
  pos_type N = cst.csa.size();
  bit_vector refbv(N, 0);

  pos_type currentIndex = cst.csa.isa[0];
  for(pos_type i = 0; i < refSize; i++) {
    refbv[currentIndex] = 1;
    currentIndex = cst.csa.psi[currentIndex];
  }

  store_to_file(refbv, refbvFilePath);
}

int main(int argc, char **argv){

  if(argc < 5){
    cerr << "Few arguments." << endl;
    cerr << "Usage : ./construct_index <ref.fa> <filtered-read.fastq> <index-file-name> <thread-num>" << endl;
    return 0;
  }

  string referenceFilePath = string(argv[1]);
  string readFilePath = string(argv[2]);
  string tmpFilePath = readFilePath + "_reference_read_text.tmp";
  ifstream ifsReference(referenceFilePath), ifsRead(readFilePath);
  ofstream ofsTmp(tmpFilePath);

  string input, referenceString;
  pos_type referenceSize = 0;
  while(getline(ifsReference, input)) {
    if(input[0] == '>') {
      referenceString += termCharacter;
      referenceSize++;
      continue;
    }
    referenceString += input;
    referenceSize += input.size();
  }
  referenceString += termCharacter;
  referenceSize++;
  ifsReference.close();

  ofsTmp << referenceString;
  referenceString.clear();

  string name, seq, blankline, qual;
  while(getline(ifsRead, name)) {
    getline(ifsRead, seq);
    getline(ifsRead, blankline);
    getline(ifsRead, qual);

    if(!isValidRead(seq)) continue;

    ofsTmp << seq << termCharacter;
    string rcseq = getRCRead(seq);
    ofsTmp << rcseq << termCharacter;
  }
  ifsRead.close();
  ofsTmp.close();

  string indexFilePath = string(argv[3]);
  string CSTFilePath = indexFilePath + ".cst";
  string refbvFilePath = indexFilePath + ".refbv";

  cerr << "[construct_index] ------ Index construction ------ " << endl;

  CST cst;

  memory_monitor::start();
  auto start = timer::now();

  construct(cst, tmpFilePath, 1);
  store_to_file(cst, CSTFilePath);

  auto stop = timer::now();
  cerr << "[construct_index] CST construction : " << duration_cast<seconds>(stop-start).count() << " seconds." << endl;

  construct_refbv(cst, referenceSize, refbvFilePath);

  stop = timer::now();
  cerr << "[construct_index] RefBitvector contruction : " << duration_cast<seconds>(stop-start).count() << " seconds." << endl;

  memory_monitor::stop();
  cerr << "[construct_index] peak usage = " << memory_monitor::peak() / (1024*1024) << " MB" << endl;

  return 0;
}
