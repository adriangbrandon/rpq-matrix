//
// Created by adrian on 8/5/25.
//

//
// Created by adrian on 8/5/25.
//
#include <iostream>
#include <sstream>
#include <chrono>
#include <vector>
#include <bm_k2_tree.hpp>
#include <file.hpp>
#include <algorithm>

#define NANO_TO_MILLI 1000000.0


/* for experiments set both constants below to 0 */
#define CHECKRESULTS 0   //if set to 1, it checks the results of the range-queries against the expected ones.
#define VERBOSETESTS 0   //if set to 1, it prints the results of each range-query launched.


namespace rangeGenerator {
    /* see also prepare-range-queries.cpp */
    /******* Reads Random Window Queries *****************/
   class reader {

    public:

       /* Gets coordinates of window and number of expected results*/
       static bool readNextWindow(std::ifstream &input, uint32_t &r1, uint32_t &r2, uint32_t &c1, uint32_t &c2, uint & nresults) {
            std::string line;
            if (std::getline(input, line)) {
                std::stringstream ss(line);
                char sep1; //#
                //line format:    r1 r2 c1 c2 #nresults
                ss >> r1 >> r2 >> c1 >> c2 >> sep1 >> nresults;
                return true;
            }
            return false;
        }

        /* Loads nqueries windows for range-queries from a file .query.res
         * Returns the max-number of results for each window (so that we can preallocate a uint64_t buffer for the results
         * of any possible range-query in the query-set.
         */
        static uint64_t read_window_queries(const std::string &file_queries,
                                    std::vector<std::vector<uint32_t>> &queries){
            std::ifstream input_query(file_queries+ ".query");  //concatenar .query  y a otro fichero .log

            std::vector<std::vector<uint32_t>> Q;
            uint32_t r1, r2, c1, c2, nresults;
            uint64_t tres=0;

           uint64_t max_nresults = 0;
            while ( readNextWindow(input_query, r1, r2, c1, c2, nresults) ) {
                std::vector<uint32_t> window;
                window.push_back(r1);window.push_back(r2);window.push_back(c1);window.push_back(c2);window.push_back(nresults);
                tres += nresults;
                if (nresults>max_nresults) max_nresults = nresults;
                Q.push_back(window);   //r1, r2, c1, c2, nresults
            }
            input_query.close();
            queries = Q;

            uint32_t nqueries = Q.size();
            std::cout << "\t *" << nqueries <<"* Window queries Read from '" << file_queries+ ".query" << "':: Avg-points=" << tres/nqueries << "."<< std::endl;
           return max_nresults;
        }



       /* Reads a line including "nresults" points, and recovers those points*/
       static std::vector<std::pair<uint32_t, uint32_t>> readNextWindowPoints(std::ifstream &input) {
           std::vector<std::pair<uint32_t, uint32_t>> points;
           uint64_t nresults;
           std::string line;

           if (std::getline(input, line)) {
               std::stringstream ss(line);
               //line format:   nresults # (r1,c1) (r2,c2) ...
               char sep1;
               ss >> nresults >> sep1;
               for (auto i=0; i<nresults; ++i ) {
                   uint32_t r,c;
                   char sepL,sepR,sepC; //  (r1,c1) (r2,c2) ...
                   ss >> sepL>> r >>sepC >> c >> sepR;
                   //points.push_back(r);
                   //points.push_back(c);
                   points.emplace_back(r, c);
               }
           }
           return points;
       }

       /* Loads the actual points (in a line) from a .query.res file which corresponds to the same line
        * from a .query file including the <r1, r2, c1, c2, nresults> describing the limits of the window query
        * */
       static void read_window_queries_PointSets(const std::string &file_queries,  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &pointSets){
            std::ifstream input_results(file_queries+ ".query.res");  //concatenar .query  y a otro fichero .log

            std::vector<std::vector<std::pair<uint32_t, uint32_t>>> QPS;
            uint64_t tres=0;
            std::vector<std::pair<uint32_t, uint32_t>> points = readNextWindowPoints(input_results);
            while ( ! points.empty() ) {
                tres += points.size();
                QPS.push_back(points);   //x, y, x, y ...
                points = readNextWindowPoints(input_results);
            }
            input_results.close();
            pointSets = QPS;

            uint32_t nqueries = QPS.size();
            std::cout << "\t *" << nqueries <<"* Window queries Read from '" << file_queries+ ".query.res" << "':: Avg-points=" << tres/nqueries << "."<< std::endl;
        }

       static void print_window_query_points(std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &pointSets, uint32_t i) {
           std::cout << "["<<pointSets[i].size() << "]:: " ;
           for (uint64_t j=0; j < pointSets[i].size(); j++) {
               std::cout << "(" << pointSets[i][j].first << "," << pointSets[i][j].second  << ")" << " ";
           }
           std::cout << std::endl;
       }

   };
}




typedef std::vector<std::vector<bool>> matrix_type;

void print_points(std::vector<std::pair<uint32_t, uint32_t>> &vp) {
    std::sort(vp.begin(), vp.end());
    for (auto & p : vp) {
        std::cout << "Point: (" << p.first << ", " << p.second << ")" << std::endl;
    }
}

void summ (matrix_type &r, matrix_type &a, matrix_type &b){
    auto n = a.size();
    auto m = a[0].size(); //    auto p = b[0].size();
    r.resize(n);
    for(auto & row : r) {
        row.resize(m);
        //row.resize(p);
    }
    for(auto i = 0; i < n; ++i) {
        for(auto j = 0; j < m; ++j) {
            r[i][j] = a[i][j] || b[i][j];
        }
    }
}

void points2matrix(std::vector<std::pair<uint32_t, uint32_t>> &points, matrix_type &m) {
    for(auto &p : points) {
        m[p.first][p.second] = true;
    }
}



int main(int argc, char **argv) {

    if (argc <3) {
        std::cerr << "Usage: " << argv[0] << " <matrices_dir> <queries_dir>" << std::endl;
        exit(1);
    }

    std::string matrices_dir = argv[1];
    std::string queries_dir = argv[2];

    uint64_t *results_buffer;  //preallocated buffer using 2*max_results of any query (known after queries are read).
    uint64_t count_results;
    uint64_t max_nresults =0;

    auto files = util::file::read_directory(matrices_dir);
    std::sort(files.begin(), files.end());

    auto queryfiles = util::file::read_directory(queries_dir);
    std::sort(queryfiles.begin(), queryfiles.end());

    std::vector<std::vector<std::vector<uint32_t>>> queriesForMatrices;
    #if CHECKRESULTS ==1
    std::vector<std::vector<std::vector<std::pair<uint32_t, uint32_t>>>> pointsetsForMatrices;
    #endif


    typedef bm_k2_tree::wrapper wrapper;
    typedef typename wrapper::matrix_type matrix;

    std::vector<matrix> m_matrices(files.size());
    uint64_t space = 0;
    std::cout << "Reading matrices..." << std::flush;


    FILE *f;
    uint nmatrices = 0;
    for(auto & file : files) {
        std::string file_path = matrices_dir + "/" + file;
        std::cout << "[" << file_path << "]" << std::endl;
        f = fopen(file_path.c_str(), "r");
        m_matrices[nmatrices] = wrapper::load(f);
        fclose(f);
        space += wrapper::space(m_matrices[nmatrices]);

        //Read range-queries for this matrix
        std::vector<std::vector<uint32_t>> queries;
        uint64_t tmpmax= rangeGenerator::reader::read_window_queries(queries_dir +  "/" + util::file::remove_extension(file), queries);
        if (tmpmax>max_nresults) max_nresults = tmpmax;
        queriesForMatrices.push_back(queries);


        #if CHECKRESULTS ==1
        //Read the results expected for the range-queries for this matrix
        std::vector<std::vector<std::pair<uint32_t, uint32_t>>> pointSets;
        rangeGenerator::reader::read_window_queries_PointSets(queries_dir + "/" + util::file::remove_extension(file), pointSets);
        pointsetsForMatrices.push_back(pointSets);
        #endif

        ++nmatrices;
    }
    std::cout << " done. [" << space << " B]" << std::endl;
    std::cout << " Max_results-range-query-buffer =  <" << max_nresults << " results>" << std::endl;

    results_buffer = (uint64_t *) malloc(2*max_nresults * sizeof(uint64_t));

    uint64_t sum = 0;
    uint64_t totNqueries = 0;

    
    std::cout << "Range Queries..." << std::flush;
    sum = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    for(uint z = 0; z < nmatrices; z++) {
        //cds::k2_tree<> &k2tree = m_matrices[z];
        matrix &k2tree = m_matrices[z];

        std::vector<std::vector<uint32_t>> &queries = queriesForMatrices[z];
#if CHECKRESULTS ==1
        std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &pointSets = pointsetsForMatrices[z];
#endif

        uint64_t i;
        //run all the range-queries over the k2-tree
        uint32_t nqueries = queries.size();
        uint64_t nresults = 0;
        for (i=0; i<nqueries; i++) {
            uint64_t nresultsExp = queries[i][4];
            count_results =0;
            //k2tree.range_query_buf(results_buffer, &count_results, 5, k2tree.dimensions, queries[i][0], queries[i][1], queries[i][2], queries[i][3],0, 0);
            count_results = wrapper::collect(k2tree,queries[i][0], queries[i][1], queries[i][2], queries[i][3],results_buffer);
            //count_results /= 2;
            nresults += count_results;

            if (nresultsExp != count_results) { std::cerr << "Error: expected " << nresultsExp << " results, but got " << count_results << std::endl; exit(1);}

            #if CHECKRESULTS ==1
            { //check results for k2_tree.range_query_buf
                std::pair<uint64_t, uint64_t>* pairs = reinterpret_cast<std::pair<uint64_t, uint64_t>*>(results_buffer);
                // Ordenar los pares (automáticamente ordena por first, luego por second)
                std::sort(pairs, pairs + count_results);

                for (uint64_t j=0; j<count_results; j++) {
                    if (pairs[j].first != pointSets[i][j].first || pairs[j].second != pointSets[i][j].second) {
                        std::cout << "\t \t Error: Query-line: " << queries[i][0] << " " << queries[i][1] << " " << queries[i][2] << " "
                                                                 << queries[i][3]<< " #nresults:" <<queries[i][4] << std::endl;
                        std::cout << "\t \t Error: " << count_results << " points recovered, expected " << pointSets[i].size() << std::endl;
                        std::cout << "\t \t Error: expected point is (" << pointSets[i][j].first << ","<< pointSets[i][j].second <<")"
                                            << " recovered point is (" << pairs[j].first << ","<< pairs[j].second <<") !!!" << std::endl;
                        exit(1);
                    }
                }
            }
            #endif


            //
            // matrix tmp;
            //tmp = wrapper::mult(m_matrices[i-1], m_matrices[i]);
            //sum += tmp->elems;
        }
        sum+= nresults;
        totNqueries += nqueries;
    }
    auto t2 =  std::chrono::high_resolution_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
    std::cerr << "Range queries: total = " << totNqueries << " queries, width = "<< queriesForMatrices[0][0][1] - queriesForMatrices[0][0][0] +1 << ", total results =" << sum<< std::endl;
    std::cerr << space << ";" << ns / (double) (NANO_TO_MILLI)                         << "   Total milliseconds" <<  std::endl;
    std::cerr << space << ";" << ns / (double) (nmatrices*NANO_TO_MILLI)               << "   total milliseconds per matrix (there are " << nmatrices << " matrices)"<<  std::endl;
    std::cerr << space << ";" << ns / (double) (totNqueries*NANO_TO_MILLI) << "   avg milliseconds per query " << std::endl;
    std::cout << " done. [" << sum << "]" << std::endl;

    #if CHECKRESULTS ==1
    std::cout << "\t TEST: done. [ range-queries ] "<<totNqueries << ", (width = " << queriesForMatrices[0][0][1] - queriesForMatrices[0][0][0]+1 << ")  TESTs SUCCEEDED - !" << std::endl;
    #endif

    free(results_buffer);
}
