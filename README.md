# RPQ-Matrix

Repository for the prototype source code of the paper Evaluating Regular Path Queries on Compressed Adjacency Matrices. 
### Queries and graph

The queries are available in `queries`:

- `queries/paths.tsv`: the queries of our benchmark that can be run on Wikidata dataset.

The dataset is available here: [dataset](https://zenodo.org/record/7254968). Specifically we use the file `wikidata.tar.gz`, 
which contains the data (*.dat*) and the dictionaries for subjects/objects (*.dat.SO*) and predicates (*.dat.P*).

### Instructions

To run our code, follow these steps:

```Bash
git clone https://github.com/adriangbrandon/rpq-matrix.git
cd rpq-matrix
mkdir build
cd build
cmake ..
make
```
After those steps, you can find the following executables in `build`:

- `create_pairs <dataset> <n_preds> <n_triples>`: given a path to a dataset, its number of predicates (*P*) and the number of
triples (*N*), this executable creates a file for each predicate storing the pairs that represent the positions of the 1-bits 
in the binary matrix.
- `k2_tree_build <pairs dir> <matrices dir> <logside>`: given a path to the pairs folder, an output folder, and 
the *log<sub>2</sub>(V)* where *V* is the number of vertices (side of the binary matrix), this executable creates
the *k<sup>2</sup>-tree* representation of the binary matrices. The matrices dir should end with `.k2-tree`.
- `baseline_build <pairs dir> <matrices dir> <side>`: given a path to the pairs folder, an output folder, and
the number of vertices (side of the binary matrix), this executable creates the *baseline* representation of the binary matrices.
  The matrices dir should end with `.baseline-64`.
- `baseline_32_build <pairs dir> <matrices dir> <side>`: identical to the previous one, but creates the *baseline* of 32-bits. The matrices dir should end with `.baseline-32`.
- The remaining executables are for querying each structure: `k2_tree_query`, `k2_tree_p_query`, `baseline_query`, and `baseline_32_query`.
Their parameters are `<dataset> <queries> <n_preds> <n_triples>`: a dataset, a file of queries, *P*, and *N*. The first two executables use
the same representation of the binary matrices, they only differ in the algorithms used to solve the queries.

### Example

From our Wikidata dataset we can obtain the following values:
- *N* = 958844164
- *P* = 5419
- *V* = 348945080

```Bash
mkdir pairs
cd pairs
../build/create_pairs ../dataset/wikidata-enumerated.dat 5419 958844164
cd ..
./build/k2_tree_build pairs/ dataset/wikidata-enumerated.dat.k2-tree 29
./build/baseline_build pairs/ dataset/wikidata-enumerated.dat.baseline-64 348945080
./build/baseline_32_build pairs/ dataset/wikidata-enumerated.dat.baseline-32 348945080

./build/baseline_query dataset/wikidata-enumerated.dat rpq-matrix/queries/paths.tsv 5419 958844164
./build/k2_tree_query dataset/wikidata-enumerated.dat rpq-matrix/queries/paths.tsv 5419 958844164
```

**IMPORTANT:** The folders of our representations shall be created in the same folder of the dataset because
they depend on the dictionaries for predicates and subjects/objects.

### Authors
- Diego Arroyuelo (diego.arroyuelo@uc.cl)
- Adrián Gómez-Brandón (adrian.gbrandon@udc.es)
- Gonzalo Navarro (gnavarro@dcc.uchile.cl)
