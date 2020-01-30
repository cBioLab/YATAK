#ifndef INCLUDED_RF_MAPPING_LIB
#define INCLUDED_RF_MAPPING_LIB

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/suffix_tree_algorithm.hpp>
#include <future>
#include <thread>
#include <mutex>

#include <genomeutil.hpp>

//// namespace
using namespace sdsl;
using namespace std;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

//// type definition
using CST = cst_sct3<>;
using RS = rank_support_v<>;
using cst_node_type = CST::node_type;
using cst_size_type = CST::size_type;
using pos_type = int64_t;
using read_num_type = int64_t;
using pos_range_type = pair<pos_type, pos_type>;
using path_type = pair<pos_range_type, pos_type>;

#define dbg(x) cerr<<#x<<"="<<x<<endl

struct MatchedInfo {
  bool dir; // F/T backward/forward
  cst_node_type v;
  pos_type len;
  cst_size_type pos;
  deque<char> edgeStr;
  MatchedInfo(){}
  MatchedInfo(bool _dir, cst_node_type _v, pos_type _len, cst_size_type _pos, deque<char> _edgeStr){
    dir = _dir;
    v = _v;
    len = _len;
    pos = _pos;
    edgeStr = _edgeStr;
  }
};

struct MappedPathInfo {
  path_type path;
  int clippingPos;
  int score;
  MappedPathInfo(){}
  MappedPathInfo(path_type _path, pos_type _clippingPos, int _score) {
    path = _path; clippingPos = _clippingPos; score = _score;
  }
};

struct MappingInfo {
  path_type pathL, pathR;
  int clippingPosL, clippingPosR;
  int score;
  bool isRev;
};

vector<MappingInfo> selectTwoDisjointRanges(const vector<MappingInfo>& infos) {
  const int rangeIntersecLenLimit = 50;
  if(infos.size() == 0) return vector<MappingInfo>();
  if(infos.size() == 1) return vector<MappingInfo>{infos[0]};
  function<pair<int, int>(const MappingInfo&)> getRange = [](const MappingInfo& a) {
    const int readLen = 150; // MAGIC NUMBER
    if(!a.isRev) {
      return pair<int, int>(a.clippingPosL, a.clippingPosR);
    } else {
      return pair<int, int>(readLen - a.clippingPosR, readLen - a.clippingPosL);
    }
  };
  pair<int, int> range0 = getRange(infos[0]);
  for(MappingInfo info : infos) {
    pair<int, int> currentRange = getRange(info);
    if(range0.first <= currentRange.first && currentRange.second <= range0.second) continue;
    if(currentRange.first <= range0.first && range0.second <= currentRange.second) continue;
    int intersecLen = min(range0.second, currentRange.second) - max(range0.first, currentRange.first);
    if(intersecLen > rangeIntersecLenLimit) continue;

    return vector<MappingInfo>{infos[0], info};
  }
  return vector<MappingInfo>{infos[0]};
}

//// constant
const string alphabet = "ACGT";
const int baseScore = 1;
const int unmatchPenalty = 1;
const int gapPenalty = 5;

class SuffixTree {
public:
  CST cst;
  bit_vector refbv;
  RS refrs;
  const char termCharacter = '$';
  const pos_type leastMatchNumber = 2;

  SuffixTree(const string& indexFilePath) {
    cerr << "[rf_mapping] Loading indexing ... : " << endl;

    load_from_file(cst, indexFilePath + ".cst");
    constructRefBv(indexFilePath + ".refbv");

    cerr << "[rf_mapping] Loading completed. " << endl;
  }

  void constructRefBv(const string& refbvFilePath) {
    load_from_file(refbv, refbvFilePath);
    util::init_support(refrs, &refbv);
  }

  pos_type calcReferencePosNum(const cst_node_type& v) {
    pos_type l = cst.lb(v), r = cst.rb(v);
    return refrs.rank(r + 1) - refrs.rank(l);
  }

  bool isMajorNode(const cst_node_type& v) {
    pos_type l = cst.lb(v), r = cst.rb(v);
    return (calcReferencePosNum(v) > 0) || (r - l + 1 >= leastMatchNumber);
  }

  bool moveLeftBoundToLeftAttempt(const MatchedInfo& crtinfo, char c, bool onMajorPath) {
    cst_node_type tmpnode = cst.wl(crtinfo.v, c);
    if(!onMajorPath) {
      if(tmpnode != cst.root()) return true;
      else return false;
    } else {
      if(tmpnode != cst.root() && isMajorNode(tmpnode)) return true;
      else return false;
    }
  }

  void moveLeftBoundToLeft(MatchedInfo& crtinfo, char c) {
    crtinfo.v = cst.wl(crtinfo.v, c);
    assert(crtinfo.v != cst.root());
    crtinfo.len++;
  }

  bool moveRightBoundToLeft(MatchedInfo& crtinfo, bool isTracing) {
    cst_node_type prev = crtinfo.v;

    crtinfo.v = cst.parent(crtinfo.v);
    crtinfo.len = cst.depth(crtinfo.v);

    if(isTracing) {
      pos_type preReferenceCount = calcReferencePosNum(prev);
      pos_type crtReferenceCount = calcReferencePosNum(crtinfo.v);
      if(preReferenceCount > 0 && preReferenceCount < crtReferenceCount) { // move from ref to ref
        return false;
      }
    }
    return true;
  }

  bool moveRightBoundToRightAttempt(const MatchedInfo& crtinfo, char c, bool onMajorPath) {
    cst_node_type tmpv = crtinfo.v;
    cst_size_type tmppos = crtinfo.pos;
    if(forward_search(cst, tmpv, crtinfo.len, c, tmppos) == 0) return false;
    if(!onMajorPath) {
      if(tmpv != cst.root()) return true;
      else return false;
    } else {
      if(tmpv != cst.root() && isMajorNode(tmpv)) return true;
      else return false;
    }
  }

  void moveRightBoundToRight(MatchedInfo& crtinfo, char c) {
    forward_search(cst, crtinfo.v, crtinfo.len, c, crtinfo.pos);
    assert(crtinfo.v != cst.root());
    crtinfo.edgeStr.push_back(c); crtinfo.len++;
    if(crtinfo.len == (pos_type) cst.depth(crtinfo.v)) crtinfo.edgeStr.clear();
  }

  bool moveLeftBoundToRight(MatchedInfo& crtinfo, bool isTracing) {

    assert(crtinfo.v != cst.root());

    cst_node_type prev = crtinfo.v;

    cst_node_type parv;
    if(crtinfo.len == (pos_type) cst.depth(crtinfo.v)) parv = crtinfo.v;
    else parv = cst.parent(crtinfo.v);

    //cerr << "parv : " << parv << endl;
    //cerr << "cst.sl(parv) : " << cst.sl(parv) << endl;

    if(parv == cst.root()) {
      assert(crtinfo.edgeStr.size() > 0);
      crtinfo.edgeStr.pop_front();
      crtinfo.v = cst.root();
    } else {
      crtinfo.v = cst.sl(parv);
    }
    crtinfo.len--;

    while(!crtinfo.edgeStr.empty() && !cst.is_leaf(crtinfo.v)) {
      cst_node_type ch = cst.child(crtinfo.v, crtinfo.edgeStr.front());
      pos_type gap = cst.depth(ch) - cst.depth(crtinfo.v);
      crtinfo.v = ch;
      if(crtinfo.edgeStr.size() < gap) break;
      for(pos_type i = 0; i < gap; i++) {
        crtinfo.edgeStr.pop_front();
      }
    }

    /*if(cst.lb(crtinfo.v) == 0 && crtinfo.len > 0) {
      cerr << cst.root() << endl;
      cerr << crtinfo.v << endl;
      cerr << crtinfo.len << endl;
      string e;
      for(char ccc : crtinfo.edgeStr) e += ccc;
      cerr << "edge str = " << e << endl;
      assert(false);
    }*/

    crtinfo.pos = cst.csa.isa[cst.csa[cst.lb(crtinfo.v)] + crtinfo.len - 1];

    if(isTracing) {
      pos_type preReferenceCount = calcReferencePosNum(prev);
      pos_type crtReferenceCount = calcReferencePosNum(crtinfo.v);
      if(preReferenceCount > 0 && preReferenceCount < crtReferenceCount) { // move from ref to ref
        return false;
      }
    }
    return true;
  }

  char calcNextCharLeft(const MatchedInfo& crtinfo, char queryCharacter) {
    MatchedInfo nxtinfo;
    if(moveLeftBoundToLeftAttempt(crtinfo, queryCharacter, true)) return queryCharacter;
    vector<pair<pos_type, char>> branch;
    for(char c : alphabet) {
      if(moveLeftBoundToLeftAttempt(crtinfo, c, true)) branch.push_back(pair<pos_type, char>(cst.rb(crtinfo.v) - cst.lb(crtinfo.v), c));
    }
    if(branch.size() == 0) return termCharacter;
    sort(branch.begin(), branch.end());
    return branch.back().second;
  }

  char calcNextCharRight(const MatchedInfo& crtinfo, char queryCharacter) {
    if(moveRightBoundToRightAttempt(crtinfo, queryCharacter, true)) return queryCharacter;
    vector<pair<pos_type, char>> branch;
    for(char c : alphabet) {
      if(moveRightBoundToRightAttempt(crtinfo, c, true)) branch.push_back(pair<pos_type, char>(cst.rb(crtinfo.v) - cst.lb(crtinfo.v), c));
    }
    if(branch.size() == 0) return termCharacter;
    sort(branch.begin(), branch.end());
    return branch.back().second;
  }

  /*LL psi(LL i, int direction) {
    char crtc = termCharacter;
    for(LL c = 0; c + 1 < 0xff; c++) {
      if(i >= C[direction][c] && i < C[direction][c + 1]) {
        crtc = (char) c;
        break;
      }
    }
    return bwt[direction].select(i - C[direction][crtc] + 1, crtc);
  }*/

   void showMatchInfoDetail(const MatchedInfo& crtinfo) { // for debug
    cerr << "-------------" << endl;
    cerr << "match len = " << crtinfo.len << endl;
    cerr << "match len = " << crtinfo.pos << endl;
    cerr << "match depth = " << cst.depth(crtinfo.v) << endl;
    cerr << "match range size = " << cst.rb(crtinfo.v) - cst.lb(crtinfo.v) + 1 << endl;
    string e;
    for(char ccc : crtinfo.edgeStr) e += ccc;
    cerr << "edge str = " << e << endl;

    if(cst.is_leaf(crtinfo.v)) {
      for(pos_type i = cst.lb(crtinfo.v) - 2; i <= min(cst.lb(crtinfo.v) + 20, cst.rb(crtinfo.v) + 2); i++) {
        pos_type idx = cst.csa[i], ccount = 0;
        while(ccount < 150) {
          cerr << cst.csa.text[idx];
          idx++;
          ccount++;
        }
        cerr << endl;
      }
      cerr << "==============" << endl;
    } else {
      for(pos_type i = cst.lb(crtinfo.v); i <= min(cst.lb(crtinfo.v) + 20, cst.rb(crtinfo.v)); i++) {
        pos_type idx = cst.csa[i], ccount = 0;
        while(ccount < 150) {
          cerr << cst.csa.text[idx];
          idx++;
          ccount++;
        }
        cerr << endl;
      }
      cerr << "==============" << endl;
    }
  }
};

#endif
