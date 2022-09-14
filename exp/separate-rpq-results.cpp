#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <string>


using namespace std;

#define TIME_OUT 60000000000
#define NQ 21

int main(int argc, char **argv)
{   

    if (argc < 3) {
        cerr << "  Usage: " << argv[0] << " <query-results-file> <query-types-file> <excluded-queries-file>" << endl;
        exit(1);
    }

    std::string line;
    
    uint64_t n_line = 0, q = 0;

    unordered_map<string, uint64_t> M;

    M["v /* c"] = 1;
    M["v * c"] = 2;
    M["v + c"] = 3;
    M["c * v"] = 4;
    M["c /* v"] = 5;
    M["v / c"] = 6;
    M["v / v"] = 7;
    M["v */* c"] = 8;
    M["v |* c"] = 9;
    M["v | v"] = 10;
    M["v ^ v"] = 11;
    M["v /* v"] = 12;
    M["v */*/*/*/* c"] = 13;
    M["v * v"] = 14;
    M["v /? c"] = 15;
    M["v + v"] = 16;
    M["v /+ c"] = 17;
    M["v ^/ c"] = 18;
    M["v || v"] = 19;
    M["v /^ v"] = 20;
    M["v | c"] = 21;

    std::ofstream fp[NQ+3];

    for (int i =  1; i <= NQ+2; i++) {
        fp[i].open(string(argv[1])+"-"+to_string(i), ios::out); 
        fp[i] << "\\begin{filecontents}{q" << i << "-" << argv[1] << ".dat}" << endl;
    }

    std::ofstream fp_all, fp_v_to_v, fp_c_to_v;
    fp_all.open(string(argv[1])+"-all");
    fp_v_to_v.open(string(argv[1])+"-all_v_to_v");
    fp_c_to_v.open(string(argv[1])+"-all_c_to_v");

    std::ifstream ifs_results(argv[1], std::ifstream::in);
    std::ifstream ifs_types(argv[2], std::ifstream::in);
    std::ifstream ifs_excluded(argv[3], std::ifstream::in);

    std::set<uint64_t> m_excluded;
    uint64_t eq;

    ifs_excluded >> eq;
    do {
        m_excluded.insert(eq);
	cout << "Excluding query " << eq << "..." << endl;
    } while (ifs_excluded >> eq);


    int nQ, nResults;
    int64_t qTime;
    char semicolon;

    unordered_map<uint64_t, string> nq_to_qtype;

    for (int i = 1; i <= 2110; i++) {
        getline(ifs_types, line);
        nq_to_qtype[i] = line;
    }


    for (int i = 1; i <= 1589; i++) {
        //getline(ifs_types, line);
	ifs_results >> nQ;
	ifs_results >> semicolon;
	ifs_results >> nResults;
	ifs_results >> semicolon;
	ifs_results >> qTime;

//	if (m_excluded.find(i) != m_excluded.end()) continue;

	if (nResults == -1 || nResults == -2 || qTime > TIME_OUT) qTime = TIME_OUT;

        fp_all << (float)qTime/1000000000 << endl;

        if (nq_to_qtype[nQ].at(0)=='v' and nq_to_qtype[nQ].at(nq_to_qtype[nQ].size()-1)=='v')
            fp_v_to_v << (float)qTime/1000000000 << endl;
	else
            fp_c_to_v << (float)qTime/1000000000 << endl;

        if (M.find(nq_to_qtype[nQ]) != M.end()) {
	    //if (i == 320) cout << line << endl;	
	    fp[M[nq_to_qtype[nQ]]] << (float)qTime/1000000000 << endl;
	}
	else {
	    if (nq_to_qtype[nQ].at(0) == 'c' or nq_to_qtype[nQ].at(nq_to_qtype[nQ].size()-1) == 'c')
	        fp[NQ+1] << (float)qTime/1000000000 << endl;
	    else 
		fp[NQ+2] << (float)qTime/1000000000 << endl;
	}

	//cout << nQ << ";" << nResults << ";" << qTime << endl;
    }


    for (int i =  1; i <= NQ+2; i++) {
        fp[i] << "\\end{filecontents}";
    }

    return 0;
}

