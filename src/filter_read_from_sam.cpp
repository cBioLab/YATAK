#include <iostream>
#include <vector>
#include <map>
#include <fstream>
using namespace std;

int base2num(char c) {
  const string alpha = "ACGT";
  for(int i = 0; i < 4; i++) {
    if(alpha[i] == c) {
      return i;
    }
  }
  return -1;
}

vector<string> splitLine(const string& S){
  int N = S.size();
  string buf;
  vector<string> result;
  for(int i = 0; i < N; i++){
    if(S[i] != '\t'){
      buf += S[i];
    }else{
      if(buf.size() > 0){
        result.push_back(buf);
        buf.clear();
      }
    }
  }
  if(buf.size() > 0){
    result.push_back(buf);
  }
  return result;
}

map<char, int> parseMappingInfo(const string& S) {
  if(S == "*") return map<char, int>();
  int N = S.size();
  int num = 0;
  map<char, int> result;
  for(int i = 0 ; i < N; i++){
    if(isdigit(S[i])){
      num *= 10;
      num += S[i] - '0';
    }else{
      result[S[i]] += num;
      num = 0;
    }
  }
  return result;
}

string getline() {

}

void worker_SAMFilter(int thread_idx, ) {
  const int mappingInfoColumn = 1, baseInfoColumn = 5;
  const int K = 20;

  while(getline(cin, input)) {
    if(input[0] == '@') { // skip tag
      filteredSAMFile << input << '\n';
      continue;
    }

    vector<string> currentLine = splitLine(input);
    int mappingInfo = stoi(currentLine[mappingInfoColumn]);

    if(mappingInfo >= 256) continue; // skip supplementary
    if((mappingInfo & 2) == 0) { // not proper pair mapped

      filteredSAMFile << input << '\n';
      continue;
    }

    map<char, int> baseInfo = parseMappingInfo(currentLine[baseInfoColumn]);
    int readlen = 0;
    for(auto it : baseInfo) {
      readlen += it.second;
    }

    if(baseInfo['M'] < readlen - K) { // not proper mapped as single read
      filteredSAMFile << input << '\n';
    }
  }
}

int main(int argc, char **argv) {

  if(argc < 2) {
    cerr << "Usage : ./filter_read_from_sam <filtered.sam>" << endl;
    return 0;
  }

  ofstream filteredSAMFile(argv[1]);
  mutex mtx;

  vector<thread> threads;
  for(int i = 0; i < threadNum; i++) {
    threads
    .emplace_back(worker_SAMFilter,
          i, ref(filteredSAMFile), ref(mtx)
    );
  }

  for(auto& t : threads) {
    t.join();
  }

  filteredSAMFile.close();

  return 0;
}
