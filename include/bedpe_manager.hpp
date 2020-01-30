#include <cassert>

using namespace std;
using pos_type = int64_t;

class VariantInfo {
  // .bedpe format
  // CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tTYPE\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2

public:

  vector<string> split(const string& str, char delim) const {
    vector<string> result;
    string buf;
    for(char c : str) {
      if(c == delim) {
        result.push_back(buf);
        buf.clear();
      } else {
        buf += c;
      }
    }
    if(!buf.empty()) {
      result.push_back(buf);
    }
    return result;
  }

  string CHROM_A, CHROM_B, ID, QUAL, STRAND_A, STRAND_B, TYPE, FILTER, FORMAT, SAMPLE1, SAMPLE2;
  pos_type START_A, END_A, START_B, END_B;
  vector<pair<string, string>> infoSet;

  VariantInfo() {
    CHROM_A = CHROM_B = ID = QUAL = STRAND_A = STRAND_B = TYPE = FILTER = FORMAT = SAMPLE1 = SAMPLE2 = "*";
    START_A = END_A = START_B = END_B = -1;
  }

  VariantInfo(const string& input) {
    vector<string> result;
    result = split(input, '\t');
    string INFO;

    CHROM_A = result[0];
    CHROM_B = result[3];
    ID = result[6];
    QUAL = result[7];
    STRAND_A = result[8];
    STRAND_B = result[9];
    TYPE = result[10];
    FILTER = result[11];
    INFO = result[12];
    FORMAT = result[13];
    SAMPLE1 = result[14];
    SAMPLE2 = result[15];

    auto strToPos = [](const string& str) {
      if(str == "*") return -1LL;
      else return stoll(str);
    };

    START_A = strToPos(result[1]);
    END_A = strToPos(result[2]);
    START_B = strToPos(result[4]);
    END_B = strToPos(result[5]);

    // fix
    if(TYPE != "ONESIDE-MAPPED" && TYPE != "UNMAPPED") {
      if(CHROM_A > CHROM_B || (CHROM_A == CHROM_B && START_A > START_B)) {
        swap(CHROM_A, CHROM_B); swap(START_A, START_B); swap(END_A, END_B); swap(STRAND_A, STRAND_B);
      }
    }

    infoSet = parseINFO(INFO);
  }

  string toString() const {
    string result = "";

    result += CHROM_A + "\t";

    if(START_A >= 0) result += to_string(START_A) + "\t";
    else result += "*\t";
    if(END_A >= 0) result += to_string(END_A) + "\t";
    else result += "*\t";

    result += CHROM_B + "\t";

    if(START_B >= 0) result += to_string(START_B) + "\t";
    else result += "*\t";
    if(END_B >= 0) result += to_string(END_B) + "\t";
    else result += "*\t";

    result += ID + "\t";
    result += QUAL + "\t";
    result += STRAND_A + "\t";
    result += STRAND_B + "\t";
    result += TYPE + "\t";
    result += FILTER + "\t";
    result += toStringInfoSet() + "\t";
    result += FORMAT + "\t";
    result += SAMPLE1 + "\t";
    result += SAMPLE2;

    return result;
  };

  bool isSame(const VariantInfo& x, bool onlyOneSide) const {
    const pos_type acceptableGap = 200;
    pos_type fromposL = START_A - acceptableGap, fromposR = END_A + acceptableGap;
    pos_type toposL = START_B - acceptableGap, toposR = END_B + acceptableGap;
    auto isOverlap = [](pos_type l0, pos_type r0, pos_type l1, pos_type r1) {
      return max(l0, l1) <= min(r0, r1);
    };

    if(!onlyOneSide) {
      return (CHROM_A == x.CHROM_A && CHROM_B == x.CHROM_B && isOverlap(fromposL, fromposR, x.START_A, x.END_A) && isOverlap(toposL, toposR, x.START_B, x.END_B));
    } else {
      return (CHROM_A == x.CHROM_A && isOverlap(fromposL, fromposR, x.START_A, x.END_A)) | (CHROM_B == x.CHROM_A && isOverlap(toposL, toposR, x.START_A, x.END_A))
      | (CHROM_A == x.CHROM_B && isOverlap(fromposL, fromposR, x.START_B, x.END_B)) | (CHROM_B == x.CHROM_B && isOverlap(toposL, toposR, x.START_B, x.END_B));
    }
  }

  pos_type length() const {
    if(TYPE == "INS") {
      string insLenStr = extractInfo("SEQLEN");
      if(insLenStr.empty()) return 0;
      else return stoll(insLenStr);
    } else {
      return START_B - END_A + 1;
    }
  }

  void addInfo(string key, string val) {
    infoSet.push_back(pair<string, string>(key, val));
  }

  string extractInfo(const string& key) const {
    for(pair<string, string> i : infoSet) {
      if(i.first == key) {
        return i.second;
      }
    }
    return ""; // not found
  }

  string toStringInfoSet() const {
    string infoStr = "";
    for(pair<string, string> p : infoSet) {
      if(!infoStr.empty()) infoStr += ";";
      if(p.second == "") infoStr += p.first;
      else infoStr += p.first + "=" + p.second;
    }
    if(infoStr.empty()) return "*";
    else return infoStr;
  }

  vector<pair<string, string>> parseINFO(const string& INFO) {
    vector<pair<string, string>> result;
    vector<string> allInfos = split(INFO, ';');
    for(string i : allInfos) {
      vector<string> keyVal = split(i, '=');
      if(keyVal.size() == 1) result.push_back(pair<string, string>(keyVal[0], ""));
      else if(keyVal.size() == 2) result.push_back(pair<string, string>(keyVal[0], keyVal[1]));
    }
    return result;
  }
};

class VariantManager {
public:
  vector<VariantInfo> vinfo;

  VariantManager(){}

  void variantInfoReader(const string& filePath) {
    ifstream ifs(filePath);

    string input;
    getline(ifs, input); // skip header
    while(getline(ifs, input)) {
      VariantInfo vi(input);
      if(vi.FILTER != "PASS") {
        vinfo.push_back(vi);
        continue;
      }
      bool preAdded = false;
      const pos_type acceptableBNDGAP = 200;
      auto isOverlap = [](pos_type l0, pos_type r0, pos_type l1, pos_type r1) {
        return max(l0, l1) <= min(r0, r1);
      };
      for(VariantInfo previ : vinfo) {
        if(previ.FILTER != "PASS") continue;
        pos_type pre_start_a = previ.START_A - acceptableBNDGAP, pre_end_a = previ.END_A + acceptableBNDGAP;
        pos_type pre_start_b = previ.START_B - acceptableBNDGAP, pre_end_b = previ.END_B + acceptableBNDGAP;
        pos_type start_a = vi.START_A, end_a = vi.END_A;
        pos_type start_b = vi.START_B, end_b = vi.END_B;
        if(previ.CHROM_A == vi.CHROM_A && isOverlap(pre_start_a, pre_end_a, start_a, end_a) && isOverlap(pre_start_b, pre_end_b, start_b, end_b)) {
          preAdded = true;
          break;
        }
      }
      if(preAdded) continue;
      vinfo.push_back(vi);
    }

    ifs.close();
  }

  VariantManager(const string& filePath) {
    variantInfoReader(filePath);
  }

  void variantInfoWriter(const string& filePath) {
    ofstream ofs(filePath, std::ios::out);

    ofs << "#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tTYPE\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\n";
    for(VariantInfo v : vinfo) {
      ofs << v.toString() << "\n";
    }

    ofs.close();
  }
};
