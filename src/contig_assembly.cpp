#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <set>

using namespace std;

using LL = int64_t;
using P = pair<LL, LL>;

#define dbg(x) cerr<<#x"="<<x<<endl

int base2num(char c) {
  const string alpha = "ACGT";
  for(int i = 0; i < 4; i++) {
    if(alpha[i] == c || tolower(alpha[i]) == c) {
      return i;
    }
  }
  if(c == '~') return 4;
  else return -1;
}

char rcBase(char base) {
  switch(base) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default : return 0;
  }
}

string makeRCread(const string& read) {
  int readlen = read.size();
  string rc;
  for(int i = 0; i < readlen; i++) {
    rc += rcBase(read[i]);
  }
  reverse(rc.begin(), rc.end());
  return rc;
}

string modifyReference(string input) {
  string res;
  for(LL i = 0; i < input.size(); i++) {
    if(islower(input[i])) res += toupper(input[i]);
    else res += input[i];
  }
  return res;
}

inline LL getReadIndex(LL pos, int readlen) {
  return pos / (LL)(readlen + 1);
}

bool fileExists(const string& str) {
   ifstream fs(str);
   return fs.is_open();
}

string getFilename(const string& filepath) {
  int len = filepath.size();
  for(int i = len - 1; i >= 0; i--) {
    if(filepath[i] == '/'){
      return filepath.substr(i + 1, len - (i + 1));
    }
  }
  return filepath;
}

bool isValidRead(const string& read) {
  int N = read.size();
  for(int i = 0; i < N; i++){
    if(base2num(read[i]) < 0) return false;
  }
  return true;
}

double getTime(clock_t startTime){
  return (double) (clock() - startTime) / CLOCKS_PER_SEC;
}

int calcPostag(string input){
  string res;
  int idx = input.size() - 1;
  while(isdigit(input[idx])){
    res += input[idx];
    idx--;
  }
  reverse(res.begin(), res.end());
  return stoi(res);
}

vector<string> splitLine(const string& S){
  int N = S.size();
  string buf;
  vector<string> result;
  for(int i = 0; i < N; i++){
    if(33 <= S[i] && S[i] <= 126){
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

vector<pair<char, int>> parseMappingInfo(const string& S){
  int N = S.size();
  int num = 0;
  vector<pair<char, int>> result;
  for(int i = 0 ; i < N; i++){
    if(isdigit(S[i])){
      num *= 10;
      num += S[i] - '0';
    }else{
      result.push_back(pair<char, int>(S[i], num));
      num = 0;
    }
  }
  return result;
}

struct contigInfo {
  string name, seq;
  LL id;

  contigInfo(string n, string s, LL i) {
    name = n;
    seq = s;
    id = i;
  }
};

struct positionInfo {
  string chr;
  LL pos, prefixlen;
  bool dir;

  positionInfo(string c, LL p, LL pl, bool d) {
    chr = c;
    pos = p;
    prefixlen = pl;
    dir = d;
  }
};

LL calcReadLen(const vector<pair<char, int>>& mapInfo, bool onRef) {
  LL len = 0;
  if(!onRef) {
    for(pair<char, int> p : mapInfo) {
      if(p.first != 'S' && p.first != 'H' && p.first != 'D') len += p.second;
    }
  } else {
    for(pair<char, int> p : mapInfo) {
      if(p.first != 'S' && p.first != 'H' && p.first != 'I') len += p.second;
    }
  }
  return len;
}

struct UnionFind{
  vector<int> v;
  UnionFind(int n) : v(n, -1) {}
  void init(){ for(int i = 0;i < (int)v.size();i++)v[i]=-1; }
  int find(int x) { return v[x] < 0 ? x : v[x] = find(v[x]); }
  bool unite(int x, int y) {
    x = find(x); y = find(y);
    if (x == y) return false;
    if (-v[x] < -v[y]) swap(x, y);
    v[x] += v[y]; v[y] = x;
    return true;
  }
  bool root(int x) { return v[x] < 0; }
  bool same(int x, int y) { return find(x) == find(y); }
  int size(int x) { return -v[find(x)]; }
};

LL calcDistForAssembly(string a, string b, bool revFlag) {
  const int unmatchPenalty = 1;
  const int gapPenalty = 3;

  if(revFlag) {
    reverse(a.begin(), a.end()); reverse(b.begin(), b.end());
  }

  int N = a.size(), M = b.size();
  vector<vector<int>> dpTable(N + 1, vector<int>(M + 1, N * M * gapPenalty));
  for(int i = 0; i <= N; i++) {
    dpTable[i][0] = i * gapPenalty;
  }
  for(int j = 0; j <= M; j++) {
    dpTable[0][j] = j * gapPenalty;
  }
  for(int i = 1; i <= N; i++) {
    for(int j = 1; j <= M; j++) {

      int score = dpTable[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? 0 : unmatchPenalty);
      if(dpTable[i][j] > score) {
        dpTable[i][j] = score;
      }
      score = dpTable[i][j - 1] + (i == N ? 0 : gapPenalty);
      if(dpTable[i][j] > score) {
        dpTable[i][j] = score;
      }
      score = dpTable[i - 1][j] + (j == M ? 0 : gapPenalty);
      if(dpTable[i][j] > score) {
        dpTable[i][j] = score;
      }
    }
  }

  return dpTable[N][M];
}

LL calcScoreForAssemblyAll(string a, string b) {
  const int baseScore = 1;
  const int gapPenalty = -3;

  int N = a.size(), M = b.size();
  vector<vector<int>> dpTable(N + 1, vector<int>(M + 1, - N * M * gapPenalty));
  for(int i = 0; i <= N; i++) {
    dpTable[i][0] = 0;
  }
  for(int j = 0; j <= M; j++) {
    dpTable[0][j] = 0;
  }
  for(int i = 1; i <= N; i++) {
    for(int j = 1; j <= M; j++) {

      int score = dpTable[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? baseScore : -baseScore);
      if(dpTable[i][j] < score) {
        dpTable[i][j] = score;
      }
      score = dpTable[i][j - 1] + (i == N ? 0 : gapPenalty);
      if(dpTable[i][j] < score) {
        dpTable[i][j] = score;
      }
      score = dpTable[i - 1][j] + (j == M ? 0 : gapPenalty);
      if(dpTable[i][j] < score) {
        dpTable[i][j] = score;
      }
    }
  }

  return dpTable[N][M];
}

int main(int argc, char **argv){

  if(argc < 3) {
    cerr << "Few arguments." << endl;
    cerr << "Usage : ./contig_assembly <contig.fastq> <contig.sam> <result.fastq>" << endl;
    return 0; 
  }

  cerr << "[contig_assembly] ----- Making contigs ----- " << endl;

  string fastqFilePath = string(argv[1]);
  string samFilePath = string(argv[2]);
  string resultFilePath = string(argv[3]);
  ifstream ifsFASTQ(fastqFilePath), ifsSAM(samFilePath);
  ofstream ofsRes(resultFilePath);

  vector<contigInfo> contigs;
  int status = 0;
  LL currentReadcount = 0;
  string input, name;
  while(getline(ifsFASTQ, input)) {
    if(status == 0){
      name = input.substr(1, input.size() - 1);
    }else if(status == 1){
      contigs.push_back(contigInfo(name, input, currentReadcount));
      currentReadcount++;
    }
    status++;
    status %= 4;
  }
  ifsFASTQ.close();

  const LL nameColumn = 0;
  const LL flagColumn = 1;
  const LL chromosomeInfoColumn = 2;
  const LL posInfoColumn = 3;
  const LL mappingInfoColumn = 5;

  map<string, vector<positionInfo>> eachContigPosInfo;
  LL finalClusterNum = 0;
  const int SVSizeThreashold = 30;
  while(getline(ifsSAM, input)) {
    if(input[0] == '@') continue;

    vector<string> info = splitLine(input);
    string contigName = info[nameColumn];
    string chromosomeInfo = info[chromosomeInfoColumn];

    if(chromosomeInfo == "*") { // output directly unmapped read
      continue;
    }

    LL flag = stoll(info[flagColumn]);
    bool dirInfo;
    if(flag & 16LL) {
      dirInfo = true;
    } else {
      dirInfo = false;
    }

    LL startPosition = stoll(info[posInfoColumn]);
    vector<pair<char, int>> mapInfo = parseMappingInfo(info[mappingInfoColumn]);

    LL readLenOnRef = calcReadLen(mapInfo, true);
    LL readLenOnContig = calcReadLen(mapInfo, false);

    char firstMappedCondition = mapInfo.front().first, lastMappedCondition = mapInfo.back().first;
    LL firstMappedConditionLen = mapInfo.front().second, lastMappedConditionLen = mapInfo.back().second;
    LL startPositionOnRef = startPosition - (firstMappedCondition == 'S' || firstMappedCondition == 'H' ? firstMappedConditionLen : 0);
    LL endPositionOnRef = startPosition + readLenOnRef + (lastMappedCondition == 'S' || lastMappedCondition == 'H' ? lastMappedConditionLen : 0) - 1;
    LL actualContigLen = readLenOnContig + (firstMappedCondition == 'S' || firstMappedCondition == 'H' ? firstMappedConditionLen : 0) + (lastMappedCondition == 'S' || lastMappedCondition == 'H' ? lastMappedConditionLen : 0);

    if((firstMappedCondition == 'S' || firstMappedCondition == 'H') && (firstMappedConditionLen >= SVSizeThreashold)) {
      eachContigPosInfo[contigName].push_back(positionInfo(chromosomeInfo, startPositionOnRef + firstMappedConditionLen, firstMappedConditionLen, dirInfo));
    }
    if((lastMappedCondition == 'S' || lastMappedCondition == 'H') && (lastMappedConditionLen >= SVSizeThreashold)) {
      eachContigPosInfo[contigName].push_back(positionInfo(chromosomeInfo, endPositionOnRef - lastMappedConditionLen + 1, actualContigLen - lastMappedConditionLen, dirInfo));
    }
    LL currentPositionOnRef = startPositionOnRef, currentContigPrefixLen = 0;
    for(pair<char, int> mi : mapInfo) {
      if((mi.first == 'D' || mi.first == 'I') && mi.second >= SVSizeThreashold) {
        eachContigPosInfo[contigName].push_back(positionInfo(chromosomeInfo, currentPositionOnRef, currentContigPrefixLen, dirInfo));
        if(mi.first == 'D') {
          eachContigPosInfo[contigName].push_back(positionInfo(chromosomeInfo, currentPositionOnRef + mi.second, currentContigPrefixLen, dirInfo));
        } else if(mi.first == 'I') {
          eachContigPosInfo[contigName].push_back(positionInfo(chromosomeInfo, currentPositionOnRef, currentContigPrefixLen + mi.second, dirInfo));
        }
      }
      if(mi.first != 'I') currentPositionOnRef += mi.second;
      if(mi.first != 'D') currentContigPrefixLen += mi.second;
    }
  }
  ifsSAM.close();

  LL contigNum = contigs.size();
  UnionFind contigCluster(contigNum);
  const LL distThreashold = 5;
  for(contigInfo contig1 : contigs) {
    for(contigInfo contig2 : contigs) {
      if(contig1.id == contig2.id) continue;

      auto isSame = [](const positionInfo& a, const positionInfo& b) {
        return (a.chr == b.chr && a.pos == b.pos);
      };

      for(positionInfo pos1 : eachContigPosInfo[contig1.name]) {
        for(positionInfo pos2 : eachContigPosInfo[contig2.name]) {

          if(isSame(pos1, pos2) && !contigCluster.same(contig1.id, contig2.id)) {

            string seq1 = contig1.seq;
            if(pos1.dir) {
              seq1 = makeRCread(seq1);
            }
            string seq2 = contig2.seq;
            if(pos2.dir) {
              seq2 = makeRCread(seq2);
            }

            string prefix1 = seq1.substr(0, pos1.prefixlen), suffix1 = seq1.substr(pos1.prefixlen, contig1.seq.size() - pos1.prefixlen);
            string prefix2 = seq2.substr(0, pos2.prefixlen), suffix2 = seq2.substr(pos2.prefixlen, contig2.seq.size() - pos2.prefixlen);

            if(calcDistForAssembly(prefix1, prefix2, true) <= distThreashold && calcDistForAssembly(suffix1, suffix2, false) <= distThreashold) {
              contigCluster.unite(contig1.id, contig2.id);
              break;
            }
          }
        }
      }
    }
  }
  
  vector<vector<contigInfo>> eachClusterInfos(contigNum);
  for(contigInfo contig : contigs) eachClusterInfos[contigCluster.find(contig.id)].push_back(contig);

  string mafftInputFilePath = fastqFilePath + "_cluster_contigs.fasta";
  string mafftOutputFilePath = fastqFilePath + "_mafft_proccessed.fasta";
  string mafftLogFilePath = fastqFilePath + "_mafftlog.txt";
  for(vector<contigInfo> currentCluster : eachClusterInfos) {
    LL currentClusterSize = currentCluster.size();
    if(currentClusterSize == 0) continue;

    if(currentClusterSize == 1) {
      ofsRes << "@contig_" << finalClusterNum << '\n';
      ofsRes << currentCluster[0].seq << '\n';
      ofsRes << "+" << '\n';
      ofsRes << string(currentCluster[0].seq.size(), 'I') << '\n';
      finalClusterNum++;
      continue;
    }

    ofstream ofsTmp(mafftInputFilePath);
    for(contigInfo contig : currentCluster) {
      string rcseq = makeRCread(contig.seq);
      if(calcScoreForAssemblyAll(currentCluster[0].seq, rcseq) > calcScoreForAssemblyAll(currentCluster[0].seq, contig.seq)) {
        ofsTmp << ">" << contig.name << '\n';
        ofsTmp << rcseq << '\n';
      } else {
        ofsTmp << ">" << contig.name << '\n';
        ofsTmp << contig.seq << '\n';
      }
    }
    ofsTmp.close();

    string mafftCommand = "mafft " + mafftInputFilePath + " > " + mafftOutputFilePath + " 2> " + mafftLogFilePath;
    system(mafftCommand.c_str());

    ifstream ifsTmp(mafftOutputFilePath);
    vector<string> alignments;
    getline(ifsTmp, input);
    while(1) {
      string seq;
      while(getline(ifsTmp, input)) {
        if(input[0] == '>') break;
        seq += input;
      }
      if(seq.size() == 0) break;
      alignments.push_back(seq);
    }

    LL constructedContigLen = alignments[0].size();
    string constructedContig;
    for(LL i = 0; i < constructedContigLen; i++) {
      char mostFreqChar = '-';
      LL mostFreqCharCount = 0;
      map<char, LL> charCount;
      for(LL j = 0; j < currentClusterSize; j++) {
        if(alignments[j][i] == '-') continue;
        charCount[alignments[j][i]]++;
        if(charCount[alignments[j][i]] > mostFreqCharCount) {
          mostFreqChar = alignments[j][i];
          mostFreqCharCount = charCount[alignments[j][i]];
        }
      }
      constructedContig += toupper(mostFreqChar);
    }

    string lastContigName = "@contig_" + to_string(finalClusterNum);
    ofsRes << lastContigName << '\n';
    ofsRes << constructedContig << '\n';
    ofsRes << "+" << '\n';
    ofsRes << string(constructedContig.size(), 'I') << '\n';

    finalClusterNum++;
  }
  ofsRes.close();
  cerr << "[contig_assembly] Assembly process completed. " << endl;
  cerr << "[contig_assembly] Final contig num : " << finalClusterNum << endl;

  return 0;
}
