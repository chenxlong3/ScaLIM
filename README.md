# ScaLIM
This repo gives the code of the algorithms of the paper "Scalable Link Recommendation for Influence Maximization".

# File Structure
```
.
|-- CMakeLists.txt
|-- README.md
|-- data
|   `-- GRQC
|-- exec_files
|   |-- ScaLIM.cpp
|   |-- evaluate.cpp
|   `-- format_graph.cpp
|-- include
|   |-- CommonFuncs.h
|   |-- Experiments.h
|   |-- GraphBase.h
|   |-- IOcontroller.h
|   |-- InfGraph.h
|   |-- ResultInfo.h
|   |-- SFMT
|   |-- Timer.h
|   |-- headers.h
|   `-- serialize.h
`-- scripts
    |-- GRQC.sh
    `-- build.sh
```


# Usage
## Cmake Compile
```
cd ./scripts
./build.sh
```
## Format Graph
In the folder of one dataset, two files should be included - `attr.txt` and `edgelist_ic.txt`.
`attr.txt` gives the number of nodes and edges in the following format.

`edgelist_ic.txt` stores the edge list of the graph. An example is shown below.
```
3466	937 0.1
3466	5233    0.1
3466	826 0.1
```

```
cd ../build
./format_graph -dataset ../data/GRQC -filename edgelist_ic.txt
```