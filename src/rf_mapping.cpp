#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <future>
#include <thread>
#include <mutex>

#include "rf_mapping_lib.hpp"

MappedPathInfo calcExtensionFromMEM(const string& query, SuffixTree& suffixtree, const pos_type& matchPosition, MatchedInfo crtInfo, const bool direction, const string& tag) { // backward extension
  const int indelSize = 4;
  const int matchLengthLimit = 30;
  const pos_type MEMMatchLengthLimit = 30;
  int querySize = query.size();

  if(!direction) { // left extension
    vector<int> levdist(indelSize * 2 + 1, 0);

    for(int i = -indelSize; i <= indelSize; i++) {
      levdist[i + indelSize] = -gapPenalty * (int) abs(i);
    }

    MatchedInfo maxScoreMappedPosInfo = crtInfo;
    int pos = matchPosition - 1, maxScore = 0, maxScoreClippingPosition = matchPosition;

    while(maxScoreMappedPosInfo.v != suffixtree.cst.root() && maxScoreMappedPosInfo.len > MEMMatchLengthLimit) {
      suffixtree.moveRightBoundToLeft(maxScoreMappedPosInfo, false);
    }

    path_type maxScorePath = path_type(pos_range_type(suffixtree.cst.lb(maxScoreMappedPosInfo.v), suffixtree.cst.rb(maxScoreMappedPosInfo.v)), maxScoreMappedPosInfo.len);

    int maxScorePointer = indelSize;

    while(pos >= 0) {

      if(pos - (maxScorePointer - indelSize) < 0 || pos - (maxScorePointer - indelSize) >= querySize) break;

      if(levdist[maxScorePointer] + (pos + 1) * baseScore < maxScore) break;

      char nextCharacter = suffixtree.calcNextCharLeft(crtInfo, query[pos - (maxScorePointer - indelSize)]);

      if(nextCharacter == suffixtree.termCharacter) {
        if(crtInfo.v == suffixtree.cst.root() || !suffixtree.moveRightBoundToLeft(crtInfo, true) || crtInfo.len < matchLengthLimit) {
          break;
        }
        suffixtree.moveRightBoundToLeft(maxScoreMappedPosInfo, false);
        continue;
      }

      suffixtree.moveLeftBoundToLeft(crtInfo, nextCharacter);
      suffixtree.moveLeftBoundToLeft(maxScoreMappedPosInfo, nextCharacter);

      vector<int> nxtlevdist(indelSize * 2 + 1, 0);
      for(int k = 0; k < indelSize * 2 + 1; k++) {
        if(pos - (k - indelSize) < 0 || pos - (k - indelSize) >= querySize) {
          nxtlevdist[k] = levdist[k];
          continue;
        }
        int val = levdist[k] + (query[pos - (k - indelSize)] != nextCharacter ? -baseScore : +baseScore);
        if(k > 0){
          val = max(val, nxtlevdist[k - 1] - gapPenalty);
        }
        if(k + 1 < indelSize * 2 + 1){
          val = max(val, levdist[k + 1] - gapPenalty);
        }
        nxtlevdist[k] = val;
      }

      levdist = nxtlevdist;
      if(*max_element(levdist.begin(), levdist.end()) > maxScore) {
        maxScore = *max_element(levdist.begin(), levdist.end());
        for(int k = 0; k < 2 * indelSize + 1; k++) {
          if(levdist[k] == maxScore) {
            maxScorePointer = k;
            maxScoreClippingPosition = pos - (k - indelSize);
            while(maxScoreMappedPosInfo.v != suffixtree.cst.root() && maxScoreMappedPosInfo.len > MEMMatchLengthLimit) {
              suffixtree.moveRightBoundToLeft(maxScoreMappedPosInfo, false);
            }
            maxScorePath = path_type(pos_range_type(suffixtree.cst.lb(maxScoreMappedPosInfo.v), suffixtree.cst.rb(maxScoreMappedPosInfo.v)), maxScoreMappedPosInfo.len);
          }
        }
      }
      pos--;
    }

    return MappedPathInfo(maxScorePath, maxScoreClippingPosition, maxScore);
  } else { // right extension

    // get corresponding reverse range
    MatchedInfo rcrtInfo = MatchedInfo(true, suffixtree.cst.root(), 0, 0, deque<char>());
    for(int i = matchPosition - crtInfo.len; i < matchPosition; i++) {
      suffixtree.moveRightBoundToRight(rcrtInfo, query[i]);
    }
    crtInfo = rcrtInfo;

    vector<int> levdist(indelSize * 2 + 1, 0);

    for(int i = -indelSize; i <= indelSize; i++) {
      levdist[i + indelSize] = -gapPenalty * (int) abs(i);
    }

    MatchedInfo maxScoreMappedPosInfo = crtInfo;
    int pos = matchPosition, maxScore = 0, maxScoreClippingPosition = matchPosition;

    while(maxScoreMappedPosInfo.v != suffixtree.cst.root() && maxScoreMappedPosInfo.len > MEMMatchLengthLimit) {
      suffixtree.moveLeftBoundToRight(maxScoreMappedPosInfo, false);
    }

    path_type maxScorePath = path_type(pos_range_type(suffixtree.cst.lb(maxScoreMappedPosInfo.v), suffixtree.cst.rb(maxScoreMappedPosInfo.v)), maxScoreMappedPosInfo.len);

    int maxScorePointer = indelSize;

    while(pos < querySize) {

      if(pos + (maxScorePointer - indelSize) < 0 || pos + (maxScorePointer - indelSize) >= querySize) break;

      if(levdist[maxScorePointer] + (querySize - pos) * baseScore < maxScore) break;

      char nextCharacter = suffixtree.calcNextCharRight(crtInfo, query[pos + (maxScorePointer - indelSize)]);

      if(nextCharacter == suffixtree.termCharacter) {
        if(crtInfo.v == suffixtree.cst.root() || !suffixtree.moveLeftBoundToRight(crtInfo, true) || crtInfo.len < matchLengthLimit) {
          break;
        }
        if(maxScoreMappedPosInfo.v != suffixtree.cst.root()) suffixtree.moveLeftBoundToRight(maxScoreMappedPosInfo, false);
        continue;
      }

      suffixtree.moveRightBoundToRight(crtInfo, nextCharacter);
      suffixtree.moveRightBoundToRight(maxScoreMappedPosInfo, nextCharacter);

      vector<int> nxtlevdist(indelSize * 2 + 1, 0);
      for(int k = 0; k < indelSize * 2 + 1; k++) {
        if(pos + (k - indelSize) < 0 || pos + (k - indelSize - k) >= querySize) {
          nxtlevdist[k] = levdist[k];
          continue;
        }
        int val = levdist[k] + (query[pos + (k - indelSize)] != nextCharacter ? -baseScore : +baseScore);
        if(k > 0){
          val = max(val, nxtlevdist[k - 1] - gapPenalty);
        }
        if(k + 1 < indelSize * 2 + 1){
          val = max(val, levdist[k + 1] - gapPenalty);
        }
        nxtlevdist[k] = val;
      }

      levdist = nxtlevdist;
      if(*max_element(levdist.begin(), levdist.end()) > maxScore) {
        maxScore = *max_element(levdist.begin(), levdist.end());
        for(int k = 0; k < 2 * indelSize + 1; k++) {
          if(levdist[k] == maxScore) {
            maxScorePointer = k;
            maxScoreClippingPosition = pos + (k - indelSize) + 1;
            while(maxScoreMappedPosInfo.v != suffixtree.cst.root() && maxScoreMappedPosInfo.len > MEMMatchLengthLimit) {
              suffixtree.moveLeftBoundToRight(maxScoreMappedPosInfo, false);
            }
            maxScorePath = path_type(pos_range_type(suffixtree.cst.lb(maxScoreMappedPosInfo.v), suffixtree.cst.rb(maxScoreMappedPosInfo.v)), maxScoreMappedPosInfo.len);
          }
        }
      }
      pos++;
    }

    return MappedPathInfo(maxScorePath, maxScoreClippingPosition, maxScore);
  }
}

bool MEMFilter(const string& query, SuffixTree& suffixtree, const string& tag) {
  int querySize = query.size();
  const int MEMMatchLengthLimit = 130;

  MatchedInfo crtInfo = (MatchedInfo){ false, suffixtree.cst.root(), 0, 0, deque<char>() };

  // Filter 1 : exact match of enough length
  for(int i = querySize - 1; i >= 0; i--) {
    MatchedInfo nxtInfo;
    while(!suffixtree.moveLeftBoundToLeftAttempt(crtInfo, query[i], false)) {
      suffixtree.moveRightBoundToLeft(crtInfo, false);
    }
    suffixtree.moveLeftBoundToLeft(crtInfo, query[i]);

    pos_type currentMatchPositonLeft = i;
    pos_type currentMatchPositonRight = i + crtInfo.len;

    if(crtInfo.len >= MEMMatchLengthLimit) {
      return false;
    }
  }

  return true;
}

vector<MappingInfo> calcBreakPosition(const string& query, SuffixTree& suffixtree, const string& tag, bool isRev, bool& foundWellMappedPath) {
  const pos_type MEMMatchLengthLimit = 30;
  const pos_type mappingScoreLimit = 130;
  int querySize = query.size();

  MatchedInfo crtInfo = MatchedInfo(false, suffixtree.cst.root(), 0, 0, deque<char>());
  vector<MatchedInfo> MEMInfo;
  vector<pos_type> MEMMatchPosition;

  // calc all of MEM
  for(int i = querySize - 1; i >= 0; i--) {
    bool savedMEM = false;
    while(!suffixtree.moveLeftBoundToLeftAttempt(crtInfo, query[i], true)) { // cannot extend to left
      if(!savedMEM) {
        savedMEM = true;
        if(crtInfo.len >= MEMMatchLengthLimit) {
          MEMInfo.push_back(crtInfo);
          MEMMatchPosition.push_back(i + 1);
        }
      }
      suffixtree.moveRightBoundToLeft(crtInfo, false);
    }
    suffixtree.moveLeftBoundToLeft(crtInfo, query[i]);
  }

  if(crtInfo.len >= MEMMatchLengthLimit) {
    MEMInfo.push_back(crtInfo);
    MEMMatchPosition.push_back(0);
  }

  if(MEMInfo.size() == 0) { // cannot find seed
    return vector<MappingInfo>();
  }

  vector<int> MEMInfoOrder(MEMInfo.size());
  iota(MEMInfoOrder.begin(), MEMInfoOrder.end(), 0);
  sort(MEMInfoOrder.begin(), MEMInfoOrder.end(), [&](const int& a, const int& b){
    return MEMInfo[a].len > MEMInfo[b].len;
  });

  vector<MappingInfo> resultInfos;
  const int resultNumLimit = 2;

  for(int currentOrder = 0; currentOrder < MEMInfo.size(); currentOrder++) {

    int crtinfoIndex = MEMInfoOrder[currentOrder];

    MatchedInfo maxInfo = MEMInfo[crtinfoIndex];
    int matchPosition = MEMMatchPosition[crtinfoIndex];

    int l = matchPosition, r = matchPosition + maxInfo.len;

    MappedPathInfo pathInfoLeft = calcExtensionFromMEM(query, suffixtree, l, maxInfo, false, tag);
    MappedPathInfo pathInfoRight = calcExtensionFromMEM(query, suffixtree, r, maxInfo, true, tag);
    int currentScore = pathInfoLeft.score + pathInfoRight.score + baseScore * maxInfo.len;

    if(currentScore >= mappingScoreLimit) {
      foundWellMappedPath = true;
      return vector<MappingInfo>();
    }

    MappingInfo newInfo;
    newInfo.score = currentScore;

    newInfo.clippingPosL = pathInfoLeft.clippingPos;
    newInfo.pathL = pathInfoLeft.path;

    newInfo.clippingPosR = pathInfoRight.clippingPos;
    newInfo.pathR = pathInfoRight.path;

    newInfo.isRev = isRev;

    resultInfos.push_back(newInfo);
  }

  sort(resultInfos.begin(), resultInfos.end(), [](const MappingInfo& a, const MappingInfo& b) {
    return a.score > b.score;
  });

  resultInfos = selectTwoDisjointRanges(resultInfos);

  return resultInfos;
}

bool fastGetLine(ifstream& ifs, string& input) {
  if((ifs >> input)) return true;
  else return false;
}

bool inputProcess(ifstream& ifs, string& tag, string& read, string& quality, mutex& mtx){
  string blankline;
  {
    lock_guard<mutex> lock(mtx);
    if(!fastGetLine(ifs, tag)) return false;
    fastGetLine(ifs, read);
    fastGetLine(ifs, blankline);
    fastGetLine(ifs, quality);
  }
  return true;
}

void worker_Read(ifstream& ifsTumor, SuffixTree& suffixtree, vector<pair<MappingInfo, string>>& mappingInfoAndRead, read_num_type& currentReadCount, mutex& mtx) {

  const pos_type clippingLengthLimit = 20;
  const int resultNumLimit = 2;
  string tag, read, quality;

  while(inputProcess(ifsTumor, tag, read, quality, mtx)) {

    if(!isValidRead(read)) {
      continue;
    }

    currentReadCount++;

    if(currentReadCount % 100000 == 0){
      cerr << "[rf_mapping] currentReadCount = " << currentReadCount << endl;
    }

    bool isSuspicious = true;
    for(int strand = 0; strand < 2; strand++, read = getRCRead(read)) {
      // MEM filter
      if(!MEMFilter(read, suffixtree, tag)) {
        isSuspicious = false;
        break;
      }
    }

    if(isSuspicious) {
      vector<MappingInfo> mappingInfosOfRead;
      bool wellMapped = false;
      for(int strand = 0; strand < 2; strand++, read = getRCRead(read)) {
        vector<MappingInfo> currentMappingInfos = calcBreakPosition(read, suffixtree, tag, (strand == 1), wellMapped);

        if(wellMapped) {
          break;
        }

        for(int i = 0; i < currentMappingInfos.size(); i++) {
          mappingInfosOfRead.push_back(currentMappingInfos[i]);
        }
      }

      if(wellMapped) continue;

      sort(mappingInfosOfRead.begin(), mappingInfosOfRead.end(), [](const MappingInfo& x, const MappingInfo& y){
        return x.score > y.score;
      });
      mappingInfosOfRead = selectTwoDisjointRanges(mappingInfosOfRead);

      for(int i = 0; i < min(resultNumLimit, (int) mappingInfosOfRead.size()); i++) {
        int readLen = read.length();
        if(mappingInfosOfRead[i].clippingPosL < clippingLengthLimit && readLen - mappingInfosOfRead[i].clippingPosR < clippingLengthLimit) {
          wellMapped = true;
          break;
        }
      }

      if(wellMapped) continue;

      for(int i = 0; i < min(resultNumLimit, (int) mappingInfosOfRead.size()); i++) {
        if(!mappingInfosOfRead[i].isRev) {
          mappingInfoAndRead.push_back(pair<MappingInfo, string>(mappingInfosOfRead[i], read));
        } else {
          string rcRead = getRCRead(read);
          mappingInfoAndRead.push_back(pair<MappingInfo, string>(mappingInfosOfRead[i], rcRead));
        }
      }
    }
  }
}

bool isValidSVContig(const string& prefix, const string& suffix) {
  map<char, int> basecount;
  const double threasholdRate = 0.9;
  for(int i = 0; i < prefix.size(); i++) {
    basecount[prefix[i]]++;
  }
  for(auto it : basecount) {
    double p = it.second;
    double q = prefix.size();
    if(p / q > threasholdRate) {
      return false;
    }
  }
  for(auto it1 : basecount) {
    for(auto it2 : basecount) {
      if(it1.first == it2.first) continue;
      double p = it1.second + it2.second;
      double q = prefix.size();
      if(p / q > threasholdRate) {
        return false;
      }
    }
  }

  basecount.clear();
  for(int i = 0; i < suffix.size(); i++) {
    basecount[suffix[i]]++;
  }
  for(auto it : basecount) {
    double p = it.second;
    double q = suffix.size();
    if(p / q > threasholdRate) {
      return false;
    }
  }
  for(auto it1 : basecount) {
    for(auto it2 : basecount) {
      if(it1.first == it2.first) continue;
      double p = it1.second + it2.second;
      double q = suffix.size();
      if(p / q > threasholdRate) {
        return false;
      }
    }
  }

  return true;
}

int main(int argc, char **argv){

  if(argc < 5) {
    cerr << "Few arguments." << endl;
    cerr << "Usage : ./rf_mapping <tumor.fastq> <index-file> <output.fastq> <thread-num>" << endl;
    return 0;
  }

  cerr << "[rf_mapping] ----- Mapping and assembling ----- " << endl;

  string tumorFilePath;
  string indexFilePath, extractedFilePath;
  string midFilePath, invalidFilePath;
  int threadNum;

  tumorFilePath = string(argv[1]);
  indexFilePath = string(argv[2]);
  extractedFilePath = string(argv[3]);
  threadNum = stoi(string(argv[4]));

  midFilePath = extractedFilePath + ".mid.fastq";
  invalidFilePath = extractedFilePath + ".invalid.fastq";

  auto startTimer = timer::now();
  memory_monitor::start();

  SuffixTree suffixtree(indexFilePath);

  cerr << "[rf_mapping] Mapping process started." << endl;

  read_num_type suspiciousReadcount = 0, clusteredReadcount = 0, currentReadCount = 0;
  vector<vector<pair<MappingInfo, string>>> mappingInfoAndRead;
  mappingInfoAndRead.resize(threadNum);

  mutex mtx;
  ifstream ifsTumor(tumorFilePath);

  vector<thread> threads;
  for(int i = 0; i < threadNum; i++) {
    threads
      .emplace_back(worker_Read,
        ref(ifsTumor), ref(suffixtree), ref(mappingInfoAndRead[i]), ref(currentReadCount), ref(mtx)
      );
  }

  for(auto& t : threads) {
    t.join();
  }

  cerr << "[rf_mapping] Mapping process completed." << endl;
  auto stopTimer = timer::now();
  memory_monitor::stop();
  cerr << "[rf_mapping] Mapping process : " << duration_cast<seconds>(stopTimer - startTimer).count() << " seconds." << endl;
  cerr << "[rf_mapping] peak usage = " << memory_monitor::peak() / (1024*1024) << " MB" << endl;
  
  // Assembling process
  memory_monitor::start();
  cerr << "[rf_mapping] Assembling process started." << endl;
  map<path_type, vector<pair<string, string>>> Cluster[2];
  read_num_type contigIndexMid = 0;

  ofstream ofsMid(midFilePath);

  const pos_type clippingLengthLimit = 20;
  for(int threadidx = 0; threadidx < threadNum; threadidx++) {
    for(read_num_type i = 0; i < (read_num_type) mappingInfoAndRead[threadidx].size(); i++) {

      MappingInfo& currentMappingInfo = mappingInfoAndRead[threadidx][i].first;
      string& currentReadSeq = mappingInfoAndRead[threadidx][i].second;
      int currentReadLen = currentReadSeq.length();
      string prefix, suffix;

      string postag;

      if(currentMappingInfo.clippingPosL >= clippingLengthLimit && currentReadLen - currentMappingInfo.clippingPosR < clippingLengthLimit) {
        prefix = currentReadSeq.substr(0, currentMappingInfo.clippingPosL);
        suffix = currentReadSeq.substr(currentMappingInfo.clippingPosL);
        reverse(prefix.begin(), prefix.end());
        Cluster[0][currentMappingInfo.pathL].push_back(pair<string, string>(prefix, suffix));

        postag += "_L" + to_string(currentMappingInfo.clippingPosL);
      }

      if(currentReadLen - currentMappingInfo.clippingPosR >= clippingLengthLimit && currentMappingInfo.clippingPosL < clippingLengthLimit) {
        prefix = currentReadSeq.substr(0, currentMappingInfo.clippingPosR);
        suffix = currentReadSeq.substr(currentMappingInfo.clippingPosR);
        reverse(prefix.begin(), prefix.end());
        Cluster[1][currentMappingInfo.pathR].push_back(pair<string, string>(prefix, suffix));

        postag += "_R" + to_string(currentMappingInfo.clippingPosR);
      }

      ofsMid << "@contig_" + to_string(contigIndexMid++) + postag << "\n";
      ofsMid << currentReadSeq << "\n";
      ofsMid << "+" << "\n";
      ofsMid << string(currentReadLen, 'I') << "\n";
    }
  }
  ofsMid.close();

  // make contig
  read_num_type suspiciousContigcount = 0, clusterCount = 0;
  ofstream ofsExtracted(extractedFilePath), ofsInvalid(invalidFilePath);
  map<string, bool> usedContig;
  vector<pair<string, string>> obtainedContigs;
  read_num_type clusterIndex = 0, invalidContigNum = 0;
  for(int dir = 0; dir < 2; dir++) {
    for(auto clusterElem : Cluster[dir]) {
      clusterCount++;
      vector<pair<string, string>>& SVReadPair = clusterElem.second;

      if(SVReadPair.size() < 4) continue;

      vector<vector<pair<string, string>>> minCluster;

      auto calcDistCheck = [](const pair<string, string>& a, const pair<string, string>& b) {
        const int hammingDistTheashold = 4;
        int d = 0;
        int l1 = min(a.first.size(), b.first.size());
        int l2 = min(a.second.size(), b.second.size());
        for(int i = 0; i < l1; i++) {
          if(a.first[i] != b.first[i]) {
            d++;
            if(d > hammingDistTheashold) return true;
          }
        }
        for(int i = 0; i < l2; i++) {
          if(a.second[i] != b.second[i]) {
            d++;
            if(d > hammingDistTheashold) return true;
          }
        }
        return false;
      };

      minCluster.push_back(vector<pair<string, string>>{SVReadPair[0]});

      for(read_num_type j = 1; j < SVReadPair.size(); j++) {
        bool done = false;
        for(read_num_type crtClusterNum = 0; crtClusterNum < minCluster.size(); crtClusterNum++) {
          bool ok = true;
          for(read_num_type k = 0; k < minCluster[crtClusterNum].size(); k++) {
            if(calcDistCheck(SVReadPair[j], minCluster[crtClusterNum][k])) {
              ok = false;
              break;
            }
          }
          if(ok) {
            done = true;
            minCluster[crtClusterNum].push_back(SVReadPair[j]);
            break;
          }
        }
        if(!done) {
          minCluster.push_back(vector<pair<string, string>>{SVReadPair[j]});
        }
      }

      for(read_num_type crtClusterNum = 0; crtClusterNum < minCluster.size(); crtClusterNum++) {

        if(minCluster[crtClusterNum].size() < 4) continue;

        vector<pair<string, string>>& clusterVector = minCluster[crtClusterNum];

        string prefix, suffix;
        int prefixlen = 0, suffixlen = 0;
        set<int> prefixLenSet, suffixLenSet;
        for(read_num_type j = 0; j < (read_num_type) clusterVector.size(); j++) {
          prefixlen = max(prefixlen, (int) clusterVector[j].first.length());
          prefixLenSet.insert((int) clusterVector[j].first.length());
          suffixlen = max(suffixlen, (int) clusterVector[j].second.length());
          suffixLenSet.insert((int) clusterVector[j].second.length());
        }

        if(prefixLenSet.size() < 4 || suffixLenSet.size() < 4) continue;

        for(int k = 0; k < prefixlen; k++) {
          vector<read_num_type> basecnt(4, 0);
          read_num_type maxcount = 0;
          char c;
          for(read_num_type j = 0; j < clusterVector.size(); j++) {
            int nc = -1;
            if(k < clusterVector[j].first.size()) nc = getNumberFromBase(clusterVector[j].first[k]);
            if(nc >= 0 && nc < 4) {
              basecnt[nc]++;
              if(basecnt[nc] > maxcount) {
                maxcount = basecnt[nc];
                c = alphabet[nc];
              }
            }
          }
          assert(maxcount > 0);
          prefix += c;
        }
        reverse(prefix.begin(), prefix.end());

        for(int k = 0; k < suffixlen; k++) {
          vector<read_num_type> basecnt(4, 0);
          read_num_type maxcount = 0;
          char c;
          for(read_num_type j = 0; j < clusterVector.size(); j++) {
            int nc = -1;
            if(k < clusterVector[j].second.size()) nc = getNumberFromBase(clusterVector[j].second[k]);
            if(nc >= 0 && nc < 4) {
              basecnt[nc]++;
              if(basecnt[nc] > maxcount) {
                maxcount = basecnt[nc];
                c = alphabet[nc];
              }
            }
          }
          assert(maxcount > 0);
          suffix += c;
        }

        string contig = prefix + suffix;
        if(!isValidSVContig(prefix, suffix)) {
          ofsInvalid << "@contig_" + to_string(invalidContigNum++) << '\n';
          ofsInvalid << contig << '\n';
          ofsInvalid << "+" << '\n';
          ofsInvalid << string(contig.size(), 'I') << '\n';
          continue;
        }

        string rccontig = getRCRead(contig);

        if(usedContig[contig]) continue;
        usedContig[contig] = usedContig[rccontig] = true;

        ofsExtracted << "@contig_" + to_string(clusterIndex++) << '\n';
        ofsExtracted << contig << '\n';
        ofsExtracted << "+" << '\n';
        ofsExtracted << string(contig.size(), 'I') << '\n';
      }
    }
  }

  ofsExtracted.close();
  ofsInvalid.close();

  cerr << "[rf_mapping] Assembling process completed. " << endl;
  stopTimer = timer::now();
  memory_monitor::stop();
  cerr << "[rf_mapping] Assembling process : " << duration_cast<seconds>(stopTimer - startTimer).count() << " seconds." << endl;
  cerr << "[rf_mapping] peak usage = " << memory_monitor::peak() / (1024*1024) << " MB" << endl;

  return 0;
}
