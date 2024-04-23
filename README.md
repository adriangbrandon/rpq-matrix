# Ring-RPQ

Repository for the prototype source code of the paper Time- and Space-Efficient Regular Path Queries on Graphs. This is just a prototype version aiming at reproducing the experiments of the paper. A final version will be added soon.

### Queries and graph

The queries are available in `data`:

- `data/wikidata/paths.tsv`: the queries of the first benchmark and they can be run on Wikidata dataset.
- `data/wikidata/paths-split-73.tsv`: the queries used in the second benchmark on Wikidata dataset to check the performance of splitting RPQs into sub-RPQs.
- `data/wikidata/paths.tsv`: the queries of the the second benchmar on YAGO2s.

The datasets are available here: [datasets](https://zenodo.org/record/7254968). For each dataset there are two configurations:

- `<name>.tar.gz`: this file contains the data (*.dat*) and the dictionaries for subjects/objects (*.dat.SO*) and predicates (*.dat.P*).
- `<name>.bfs.tar.gz`: as the previous file, but the subject/object identifiers are ordering by traversing the graph following a BFS search.

Once any `<name>.tar.gz` or `<name>.bfs.tar.gz` is decompressed, their files have to be kept in the same path.
### Instructions

To run our code, please install an extended version of the library SDSL. Go to this [this repository](https://github.com/adriangbrandon/sdsl-lite) and follow the instructions.

After the extended version of SDSL is installed, clone this repository and follow these steps:

```Bash
git clone https://github.com/adriangbrandon/Ring-RPQ.git
cd Ring-RPQ
mkdir build
cd build
cmake ..
make

mkdir pairs
cd pairs
../bin/create_pairs ../dataset/wikidata-enumerated.dat 5419 958844164
cd ..
./bin/k2_tree_build pairs/ k2-tree/wikidata-enumerated.dat.matrices 29
./bin/baseline_build pairs/ baseline-64/wikidata-enumerated.dat.matrices 348945080
./bin/baseline_32_build pairs/ baseline-64/wikidata-enumerated.dat.matrices 348945080
```

This shall create several executables:

- `build-index-basic <path_to_dat>`: given a path to the file *.dat*, it builds in the same folder our **Ring** or **RingA** system depending on
the kind of dataset you are choosing. In case of using `<name>.tar.gz` you will obtain **Ring**. Otherwise, using the dataset of `<name>.bfs.tar.gz`, it builds the **RingA**.
- `build-index-split <path_to_dat>`: given a path to the *.dat* from a dataset `<name>.bfs.tar.gz` it builds the **RingAB** in the same folder.

**IMPORTANT:** Keep all the files in the same directory.

- `query-index-basic <path_to_dat> <path_to_queries>`: given a path to the *.dat* from a dataset and a path to the queries, it runs those queries on **Ring** or **RingA** and computes their elapsed time.
- `query-index-split <path_to_dat> <path_to_queries>`: given a path to the *.dat* from a dataset and a path to the queries, it runs those queries on **RingAB** and computes their elapsed time.