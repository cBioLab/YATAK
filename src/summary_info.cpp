#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <unordered_map>
#include <iomanip>
#include <math.h>

#include <bedpe_manager.hpp>

using namespace std;

using pos_type = int64_t;

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

vector<pair<char, pos_type>> parseMappingConditionInfo(const string& S){
  pos_type N = S.size();
  pos_type num = 0;
  vector<pair<char, pos_type>> result;
  for(pos_type i = 0 ; i < N; i++){
    if(isdigit(S[i])){
      num *= 10;
      num += S[i] - '0';
    }else{
      result.push_back(pair<char, pos_type>(S[i], num));
      num = 0;
    }
  }
  return result;
}

struct variantInfo {
  string chr1, chr2;
  pos_type startpos1, endpos1, startpos2, endpos2;
  bool dir1, dir2;
  string type;
  string seq;
  int score1, score2;

  variantInfo(string c1, pos_type sp1, pos_type ep1, bool d1, string c2, pos_type sp2, pos_type ep2, bool d2, string t, string s, int sc1, int sc2) {
    chr1 = c1; chr2 = c2;
    startpos1 = sp1; startpos2 = sp2;
    endpos1 = ep1; endpos2 = ep2;
    dir1 = d1; dir2 = d2;
    type = t;
    seq = s;
    score1 = sc1;
    score2 = sc2;
  }
};

struct mappedInfo {
  string chr;
  pos_type startPosOnRef, endPosOnRef;
  pos_type startPosOnContig, endPosOnContig;
  int flag, mapq;
  bool dir;
  vector<pair<char, pos_type>> mapConditionInfo;
};

pos_type calcReadLen(const vector<pair<char, pos_type>>& mapInfo, bool onRef) {
  pos_type len = 0;
  if(!onRef) {
    for(pair<char, pos_type> p : mapInfo) {
      if(p.first != 'S' && p.first != 'H' && p.first != 'D') len += p.second;
    }
  } else {
    for(pair<char, pos_type> p : mapInfo) {
      if(p.first != 'S' && p.first != 'H' && p.first != 'I') len += p.second;
    }
  }
  return len;
}

int main(int argc, char **argv){

  if(argc < 4) {
    cerr << "Few arguments." << endl;
    cerr << "Usage : ./summary_info <contig.fastq> <contig.sam> <result.bedpe>" << endl;
    return 0; 
  }

  cerr << "[summary_info] ------ Making summary ------ " << endl;

  string fastqFilePath = string(argv[1]);
  string samFilePath = string(argv[2]);
  string resultFilePath = string(argv[3]);
  ifstream ifsFASTQ(fastqFilePath), ifsSAM(samFilePath);

  map<string, string> contigSeq;
  map<string, bool> isUsedContig;
  int status = 0;
  string input, name;
  while(getline(ifsFASTQ, input)) {
    if(status == 0){
      name = input.substr(1);
    }else if(status == 1){
      contigSeq[name] = input;
    }
    status++;
    status %= 4;
  }
  ifsFASTQ.close();

  const int nameColumn = 0;
  const int flagColumn = 1;
  const int chromosomeInfoColumn = 2;
  const int posInfoColumn = 3;
  const int mapqInfoColumn = 4;
  const int mappingConditionInfoColumn = 5;

  map<string, vector<mappedInfo>> eachContigMappedInfo;
  while(getline(ifsSAM, input)) {
    if(input[0] == '@') continue;

    vector<string> splitInput = splitLine(input);

    mappedInfo currentInfo;
    currentInfo.flag = stoi(splitInput[flagColumn]);
    currentInfo.chr = splitInput[chromosomeInfoColumn];

    if(currentInfo.flag == 4) {
      eachContigMappedInfo[splitInput[nameColumn]].push_back(currentInfo);
      continue;
    }

    if((currentInfo.flag & 16) == 0) {
      currentInfo.dir = false;
    } else {
      currentInfo.dir = true;
    }

    currentInfo.mapq = stoi(splitInput[mapqInfoColumn]);

    currentInfo.mapConditionInfo = parseMappingConditionInfo(splitInput[mappingConditionInfoColumn]);
    currentInfo.startPosOnRef = stoll(splitInput[posInfoColumn]);

    pos_type readLenOnRef = calcReadLen(currentInfo.mapConditionInfo, true);
    currentInfo.endPosOnRef = currentInfo.startPosOnRef + readLenOnRef - 1;

    pos_type readLenOnContig = calcReadLen(currentInfo.mapConditionInfo, false);

    char firstMappedCondition = currentInfo.mapConditionInfo.front().first;
    pos_type firstMappedConditionLen = currentInfo.mapConditionInfo.front().second;
    currentInfo.startPosOnContig = (firstMappedCondition == 'S' || firstMappedCondition == 'H' ? firstMappedConditionLen : 0);
    currentInfo.endPosOnContig = currentInfo.startPosOnContig + readLenOnContig - 1;

    if(currentInfo.dir) {
      pos_type actualReadLen = 0;
      for(pair<char, pos_type> p : currentInfo.mapConditionInfo) if(p.first != 'D') actualReadLen += p.second;
      pos_type actualStartPosOnContig = actualReadLen - 1 - currentInfo.endPosOnContig;
      pos_type actualEndPosOnContig = actualReadLen - 1 - currentInfo.startPosOnContig;

      currentInfo.startPosOnContig = actualStartPosOnContig;
      currentInfo.endPosOnContig = actualEndPosOnContig;
    }

    eachContigMappedInfo[splitInput[nameColumn]].push_back(currentInfo);
  }
  ifsSAM.close();

  vector<variantInfo> vinfo;
  const int SVSizeThreashold = 30;
  const int mapConditionLenThreashold = 30;
  for(auto it : eachContigMappedInfo) {
    string contigName = it.first;
    vector<mappedInfo> infos = it.second;

    isUsedContig[contigName] = true;

    if(infos.size() == 1 && infos[0].flag == 4) { // unmapped
      vinfo.push_back(variantInfo("*", -1, -1, false, "*", -1, -1, false, "unmapped-contig", contigSeq[contigName], 0, 0));
      continue;
    }

    vector<mappedInfo> disjointInfos;
    for(int i = 0; i < infos.size(); i++) {
      bool included = false;
      for(int j = 0; j < disjointInfos.size(); j++) {
        if(disjointInfos[j].startPosOnContig <= infos[i].startPosOnContig && infos[i].endPosOnContig <= disjointInfos[j].endPosOnContig) { // included
          included = true;
          break;
        }
      }
      if(!included) {
        disjointInfos.push_back(infos[i]);
      }
    }
    infos = disjointInfos;

    if(infos.size() == 1) {
      mappedInfo info = infos[0];
      int flag = info.flag;

      string chromosomeInfo = info.chr;
      vector<pair<char, pos_type>> mapInfo = info.mapConditionInfo;

      // oneside-mapped
      if((mapInfo.front().first == 'S' || mapInfo.front().first == 'H') && mapInfo.front().second >= SVSizeThreashold) {
        vinfo.push_back(variantInfo(info.chr, info.startPosOnRef, info.endPosOnRef, false, "*", -1, -1, false, "oneside-mapped-contig", contigSeq[contigName], info.mapq, 0));
      }
      if((mapInfo.back().first == 'S' || mapInfo.back().first == 'H') && mapInfo.back().second >= SVSizeThreashold) {
        vinfo.push_back(variantInfo(info.chr, info.startPosOnRef, info.endPosOnRef, true, "*", -1, -1, false, "oneside-mapped-contig", contigSeq[contigName], info.mapq, 0));
      }

      pos_type currentPos = info.startPosOnRef;

      // small indel
      for(pair<char, pos_type> p : mapInfo) {
        if(p.first == 'D' && p.second >= SVSizeThreashold) {
          pos_type anotherStartPos = currentPos + p.second;
          vinfo.push_back(variantInfo(info.chr, info.startPosOnRef, currentPos, true, info.chr, anotherStartPos, info.endPosOnRef, false, "small-deletion", "*", info.mapq, info.mapq));
        }
        if(p.first == 'I' && p.second >= SVSizeThreashold) {
          pos_type anotherStartPos = currentPos + 1;
          vinfo.push_back(variantInfo(info.chr, info.startPosOnRef, currentPos, true, info.chr, anotherStartPos, info.endPosOnRef, false, "small-insertion", "*", info.mapq, info.mapq));
        }
        if(p.first != 'I' && p.first != 'S' && p.first != 'H') currentPos += p.second;
      }
    } else if(infos.size() >= 2) {

      sort(infos.begin(), infos.end(), [](const mappedInfo& x, const mappedInfo& y) {
        return x.startPosOnContig < y.startPosOnContig;
      });

      bool tooLargeOverlapOrGap = false;
      for(int i = 0; i < infos.size() - 1; i++){
        int overlapLen = min(infos[i].endPosOnContig, infos[i + 1].endPosOnContig) - max(infos[i].startPosOnContig, infos[i + 1].startPosOnContig);
        int gapLen = infos[i + 1].startPosOnContig - infos[i].endPosOnContig;
        if(overlapLen > mapConditionLenThreashold || gapLen > mapConditionLenThreashold) {
          tooLargeOverlapOrGap = true;
        }
      }
      if(tooLargeOverlapOrGap) {
        vinfo.push_back(variantInfo("*", -1, -1, false, "*", -1, -1, false, "unsolved-contig", contigSeq[contigName], 0, 0));
        continue;
      }

      for(int i = 0; i < infos.size() - 1; i++) {
        mappedInfo info1 = infos[i];
        mappedInfo info2 = infos[i + 1];

        pos_type pos1, pos2;
        bool dir1, dir2;
        if(!info1.dir) {
          dir1 = true;
          pos1 = info1.endPosOnRef;
        } else {
          dir1 = false;
          pos1 = info1.startPosOnRef;
        }
        if(!info2.dir) {
          dir2 = false;
          pos2 = info2.startPosOnRef;
        } else {
          dir2 = true;
          pos2 = info2.endPosOnRef;
        }

        bool swapped = false;
        if((info1.chr == info2.chr && pos1 > pos2) || info1.chr > info2.chr) {
          swap(info1, info2); swap(dir1, dir2); swap(pos1, pos2);
          swapped = true;
        }

        if(disjointInfos.size() == 2) {
          if(info1.chr == info2.chr) {
            if((!swapped && !info1.dir && !info2.dir) || (swapped && info1.dir && info2.dir)) {
              vinfo.push_back(variantInfo(info1.chr, info1.startPosOnRef, info1.endPosOnRef, dir1, info2.chr, info2.startPosOnRef, info2.endPosOnRef, dir2, "deletion", "*", info1.mapq, info2.mapq));
            } else if((!swapped && !info1.dir && info2.dir) || (!swapped && info1.dir && !info2.dir) || (swapped && info1.dir && !info2.dir) || (swapped && !info1.dir && info2.dir)) {
              vinfo.push_back(variantInfo(info1.chr, info1.startPosOnRef, info1.endPosOnRef, dir1, info2.chr, info2.startPosOnRef, info2.endPosOnRef, dir2, "inversion", "*", info1.mapq, info2.mapq));
            } else if((!swapped && info1.dir && info2.dir) || (swapped && !info1.dir && !info2.dir)) {
              vinfo.push_back(variantInfo(info1.chr, info1.startPosOnRef, info1.endPosOnRef, dir1, info2.chr, info2.startPosOnRef, info2.endPosOnRef, dir2, "duplication", "*", info1.mapq, info2.mapq));
            } else {
              vinfo.push_back(variantInfo(info1.chr, info1.startPosOnRef, info1.endPosOnRef, dir1, info2.chr, info2.startPosOnRef, info2.endPosOnRef, dir2, "breakend", "*", info1.mapq, info2.mapq));
            }
          } else {
            vinfo.push_back(variantInfo(info1.chr, info1.startPosOnRef, info1.endPosOnRef, dir1, info2.chr, info2.startPosOnRef, info2.endPosOnRef, dir2, "breakend", "*", info1.mapq, info2.mapq));
          }
        } else {
          vinfo.push_back(variantInfo(info1.chr, info1.startPosOnRef, info1.endPosOnRef, dir1, info2.chr, info2.startPosOnRef, info2.endPosOnRef, dir2, "complex", "*", info1.mapq, info2.mapq));
        }
      }
    }
  }

  // solve large insertion
  vector<variantInfo> insertSolvedInfo;
  vector<bool> usedInsertion(vinfo.size(), false);
  const pos_type insertionGapThreashold = 30;
  for(int i = 0; i < vinfo.size(); i++) {
    if(usedInsertion[i]) continue;
    if(vinfo[i].type != "oneside-mapped-contig") continue;
    for(int j = i + 1; j < vinfo.size(); j++) {
      if(usedInsertion[j]) continue;
      if(vinfo[j].type != "oneside-mapped-contig") continue;

      variantInfo info1 = vinfo[i];
      variantInfo info2 = vinfo[j];

      pos_type pos1, pos2;
      if(!info1.dir1) {
        pos1 = info1.startpos1;
      } else {
        pos1 = info1.endpos1;
      }

      if(!info2.dir1) {
        pos2 = info2.startpos1;
      } else {
        pos2 = info2.endpos1;
      }

      if(info1.chr1 == info2.chr1 && info1.dir1 != info2.dir1 && (pos_type) abs(pos1 - pos2) <= insertionGapThreashold) {
        usedInsertion[i] = usedInsertion[j] = true;
        if(!info1.dir1) {
          swap(info1, info2);
        }
        insertSolvedInfo.push_back(variantInfo(info1.chr1, info1.startpos1, info1.endpos1, info1.dir1, info2.chr1, info1.endpos1 + 1, info2.endpos1, info2.dir1, "insertion", "*", info1.score1, info2.score1));
      }
    }
  }

  for(int i = 0; i < vinfo.size(); i++) {
    if(usedInsertion[i]) continue;
    insertSolvedInfo.push_back(vinfo[i]);
  }
  vinfo = insertSolvedInfo;

  // solve dup inversion
  vector<variantInfo> inversionSolvedInfo;
  vector<bool> dupInversion(vinfo.size(), false);
  const pos_type inversionGapThreashold = 30;
  for(int i = 0; i < vinfo.size(); i++) {
    if(dupInversion[i]) continue;
    if(vinfo[i].type != "inversion") continue;
    for(int j = i + 1; j < vinfo.size(); j++) {
      if(dupInversion[j]) continue;
      if(vinfo[j].type != "inversion") continue;

      variantInfo info1 = vinfo[i];
      variantInfo info2 = vinfo[j];

      pos_type pos1a, pos2a, pos1b, pos2b;
      if(!info1.dir1) {
        pos1a = info1.startpos1;
      } else {
        pos1a = info1.endpos1;
      }
      if(!info1.dir2) {
        pos2a = info1.startpos2;
      } else {
        pos2a = info1.endpos2;
      }

      if(!info2.dir1) {
        pos1b = info2.startpos1;
      } else {
        pos1b = info2.endpos1;
      }
      if(!info2.dir2) {
        pos2b = info2.startpos2;
      } else {
        pos2b = info2.endpos2;
      }

      if(info1.dir1 != info2.dir1 && info1.dir2 != info2.dir2 && (pos_type) abs(pos1a - pos1b) <= inversionGapThreashold && (pos_type) abs(pos2a - pos2b) <= inversionGapThreashold) {
        dupInversion[j] = true;
      }
    }
  }

  for(int i = 0; i < vinfo.size(); i++) {
    if(dupInversion[i]) continue;
    inversionSolvedInfo.push_back(vinfo[i]);
  }
  vinfo = inversionSolvedInfo;

  for(auto c : contigSeq) {
    if(isUsedContig[c.first]) continue;
    vinfo.push_back(variantInfo("*", -1, -1, false, "*", -1, -1, false, "unmapped-contig", c.second, -1, -1));
  }

  sort(vinfo.begin(), vinfo.end(), [](const variantInfo& a, const variantInfo& b) {
    if(a.chr1 == "*") return false;
    if(b.chr1 == "*") return true;

    if(a.chr2 == "*" && b.chr2 != "*") return false;
    if(a.chr2 != "*" && b.chr2 == "*") return true;

    pos_type apos1 = (!a.dir1 ? a.startpos1 : a.endpos1), apos2 = (!a.dir2 ? a.startpos2 : a.endpos2);
    pos_type bpos1 = (!b.dir1 ? b.startpos1 : b.endpos1), bpos2 = (!b.dir2 ? b.startpos2 : b.endpos2);

    return make_tuple(a.chr1, apos1, a.chr2, apos2) < make_tuple(b.chr1, bpos1, b.chr2, bpos2);
  });

  // output final result
  VariantManager finalResult;
  for(variantInfo i : vinfo) {
    VariantInfo vi;

    if(i.type == "deletion" || i.type == "small-deletion") {
      vi.TYPE = "DEL";
      if(i.score1 == 0 || i.score2 == 0) {
        vi.FILTER = "REJECT";
      } else {
        vi.FILTER = "PASS";
      }

      vi.CHROM_A = vi.CHROM_B = i.chr1;
      vi.STRAND_A = "+"; vi.STRAND_B = "-";
      vi.START_A = vi.END_A = i.endpos1;
      vi.START_B = vi.END_B = i.startpos2;
      vi.addInfo("SCORE1", to_string(i.score1));
      vi.addInfo("SCORE2", to_string(i.score2));
    } else if(i.type == "insertion" || i.type == "small-insertion") {
      vi.TYPE = "INS";
      if(i.score1 == 0 || i.score2 == 0) {
        vi.FILTER = "REJECT";
      } else {
        vi.FILTER = "PASS";
      }

      vi.CHROM_A = vi.CHROM_B = i.chr1;
      vi.STRAND_A = "+"; vi.STRAND_B = "-";
      vi.START_A = vi.END_A = i.endpos1;
      vi.START_B = vi.END_B = i.startpos2;
      vi.addInfo("SCORE1", to_string(i.score1));
      vi.addInfo("SCORE2", to_string(i.score2));
    } else if(i.type == "inversion") {

      vi.TYPE = "INV";
      if(i.score1 == 0 || i.score2 == 0) {
        vi.FILTER = "REJECT";
      } else {
        vi.FILTER = "PASS";
      }

      vi.CHROM_A = vi.CHROM_B = i.chr1;
      vi.STRAND_A = "+"; vi.STRAND_B = "-";
      if(!i.dir1) {
        vi.START_A = vi.END_A = i.startpos1;
        vi.START_B = vi.END_B = i.startpos2;
      } else {
        vi.START_A = vi.END_A = i.endpos1;
        vi.START_B = vi.END_B = i.endpos2;
      }

      const pos_type smallINVLengthLimit = 50;
      if (vi.START_B - vi.START_A <= smallINVLengthLimit) { // discard small inversion contig
        continue;
      }

      vi.addInfo("SCORE1", to_string(i.score1));
      vi.addInfo("SCORE2", to_string(i.score2));
    } else if(i.type == "duplication") {
      vi.TYPE = "DUP";
      if(i.score1 == 0 || i.score2 == 0) {
        vi.FILTER = "REJECT";
      } else {
        vi.FILTER = "PASS";
      }

      vi.CHROM_A = vi.CHROM_B = i.chr1;
      vi.STRAND_A = "+"; vi.STRAND_B = "-";
      vi.START_A = vi.END_A = i.startpos1;
      vi.START_B = vi.END_B = i.endpos2;

      vi.addInfo("SCORE1", to_string(i.score1));
      vi.addInfo("SCORE2", to_string(i.score2));
    } else if(i.type == "breakend") {
      vi.TYPE = "BND";
      if(i.score1 == 0 || i.score2 == 0) {
        vi.FILTER = "REJECT";
      } else {
        vi.FILTER = "PASS";
      }

      vi.CHROM_A = i.chr1;
      if(!i.dir1) {
        vi.START_A = vi.END_A = i.startpos1;
        vi.STRAND_A = "-";
      } else {
        vi.START_A = vi.END_A = i.endpos1;
        vi.STRAND_A = "+";
      }
      vi.CHROM_B = i.chr2;
      if(!i.dir2) {
        vi.START_B = vi.END_B = i.startpos2;
        vi.STRAND_B = "-";
      } else {
        vi.START_B = vi.END_B = i.endpos2;
        vi.STRAND_B = "+";
      }
      vi.addInfo("SCORE1", to_string(i.score1));
      vi.addInfo("SCORE2", to_string(i.score2));
    } else if(i.type == "complex") {
      vi.TYPE = "BNDCMP";
      if(i.score1 == 0 || i.score2 == 0) {
        vi.FILTER = "REJECT";
      } else {
        vi.FILTER = "PASS";
      }

      vi.CHROM_A = i.chr1;
      if(!i.dir1) {
        vi.START_A = vi.END_A = i.startpos1;
        vi.STRAND_A = "-";
      } else {
        vi.START_A = vi.END_A = i.endpos1;
        vi.STRAND_A = "+";
      }
      vi.CHROM_B = i.chr2;
      if(!i.dir2) {
        vi.START_B = vi.END_B = i.startpos2;
        vi.STRAND_B = "-";
      } else {
        vi.START_B = vi.END_B = i.endpos2;
        vi.STRAND_B = "+";
      }
      vi.addInfo("SCORE1", to_string(i.score1));
      vi.addInfo("SCORE2", to_string(i.score2));
    } else if(i.type == "oneside-mapped-contig") {
      vi.TYPE = "ONESIDE-MAPPED";
      vi.FILTER = "REJECT";

      vi.CHROM_A = i.chr1;
      if(!i.dir1) {
        vi.START_A = vi.END_A = i.startpos1;
        vi.STRAND_A = "-";
      } else {
        vi.START_A = vi.END_A = i.endpos1;
        vi.STRAND_A = "+";
      }

      vi.addInfo("SEQ", i.seq);
      vi.addInfo("SCORE1", to_string(i.score1));
    } else if(i.type == "unmapped-contig" || i.type == "unsolved-contig") {
      vi.TYPE = "UNMAPPED";
      vi.FILTER = "REJECT";

      vi.addInfo("SEQ", i.seq);
    }
    if(i.score1 > 0 && i.score2 > 0) {
      vi.addInfo("CONDITION", "BOTH");
    } else if(i.score1 > 0 || i.score2 > 0) {
      vi.addInfo("CONDITION", "ONESIDE");
    } else {
      vi.addInfo("CONDITION", "NEITHER");
    }
    finalResult.vinfo.push_back(vi);
  }
  finalResult.variantInfoWriter(resultFilePath);

  cerr << "[summary_info] Making summary completed. " << endl;

  return 0;
}
