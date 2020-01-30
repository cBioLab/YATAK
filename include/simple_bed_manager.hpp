#include <functional>

using namespace std;
using pos_type = int64_t;

class SimpleVariantInfo {
  // simple bed format for HACk
  // CHROM\tSTART\tEND\tTYPE\tSEQ\tTEMPSEQLEN

private:
  vector<string> split(const string& str, char delim) {
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

public:
  string CHROM, TYPE, SEQ, TEMPSEQLEN;
  pos_type START, END;

  SimpleVariantInfo() {
    SEQ = "*";
    TEMPSEQLEN = "0";
  }

  SimpleVariantInfo(string chrom, pos_type start, pos_type end, string type, string seq, string tempseqlen) {
    CHROM = chrom;
    START = start;
    END = end;
    TYPE = type;
    SEQ = seq;
    TEMPSEQLEN = tempseqlen;
  }

  SimpleVariantInfo(const string& input) {
    vector<string> result;
    result = split(input, '\t');

    CHROM = result[0];
    TYPE = result[3];
    SEQ = result[4];
    TEMPSEQLEN = result[5];

    START = stoll(result[1]);
    END = stoll(result[2]);
  }

  string toString() const {
    string result = "";

    result += CHROM + "\t";
    result += to_string(START) + "\t";
    result += to_string(END) + "\t";
    result += TYPE + "\t";
    result += SEQ + "\t";
    result += TEMPSEQLEN + "\n";

    return result;
  };

  pos_type length() const {
    if(TYPE == "insertion") {
      return SEQ.size();
    } else {
      return END - START + 1;
    }
  }

  bool operator<(const SimpleVariantInfo& x) const {
    if(CHROM == x.CHROM) {
      return pair<pos_type, pos_type>(START, END) < pair<pos_type, pos_type>(x.START, x.END);
    } else {
      return CHROM < x.CHROM;
    }
  }
};

class SimpleVariantManager {
public:
  vector<SimpleVariantInfo> vinfo;

  SimpleVariantManager(){}

  void simpleVariantInfoReader(const string& filePath) {
    ifstream ifs(filePath);

    string input;
    while(getline(ifs, input)) {
      SimpleVariantInfo svi(input);
      vinfo.push_back(SimpleVariantInfo(input));
    }

    ifs.close();
  }

  SimpleVariantManager(const string& filePath) {
    simpleVariantInfoReader(filePath);
  }

  void simpleVariantInfoWriter(const string& filePath) {
    ofstream ofs(filePath);

    for(SimpleVariantInfo v : vinfo) {
      ofs << v.toString();
    }

    ofs.close();
  }

  void filterInfo(function<bool(const SimpleVariantInfo&)> filterFunc) {
    auto filteredIter = remove_if(vinfo.begin(), vinfo.end(), filterFunc);
    vinfo.erase(filteredIter, vinfo.end());
  }
};
