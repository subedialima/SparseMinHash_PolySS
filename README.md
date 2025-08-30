# Approximate Nearest Neighbor with Real World Geometry

This repository contains the source code for SparseMinHash-PolySS. 

This research hopes to use locality sensitive hashing to find the nearest neighbors of polygons using Jaccard Distance as our metric. 
The project is primarily written in C++ using the reentrant C API for the [GEOS library](https://libgeos.org/) for the hashing algorithms.
Data visualization is done with Python using the [Shapely](https://shapely.readthedocs.io/en/stable/manual.html) along with matplotlib.

The project currently has multiple directories.

1. `data` Where the input data is stored and the output data is written. Due to the large nature of the datasets, data is not available in this repository but can be found at the [UCR Star](https://star.cs.ucr.edu/) website. 
2. `src` The actual source code for the nearest neighbor search code, written in C++20. 


Within the code itself, each function, class, type, etc in a header file has a description of its respective parameters, return values, members, etc, which I hope will make understanding the code easier. 

# 1.1 Additional Experimental evaluation of LSH-PolySS

We are using performance metrics such as **Recall@k** and **runtime** taken by LSH-PolySS to query input polygons compared to the ground truth found by brute-force (BF). 

- **Recall@k** is the proportion of top-k approximate nearest-neighbors returned by the LSH search, which are part of the set of actual nearest top-k neighbors returned by the BF search.  
- Possible values range from 0 to 1, with higher values indicating better performance.  

In all implementation methods discussed, we will show recall for different top-k neighbors when grid cell threshold is greater than 0. Since we are refining (for verification) the candidates for each query, the time to select the top-k neighbors is roughly equal. Thus, we are reporting query time for only the top 10 neighbors.

---

## 1.1.1 Sketch Similarity of Encoded Vectors without Minhashing

In order to see the benefits of minhashing and LSH, we first compare the quality of the results using only the encoded sketch vectors. We present the recall and runtime for the sketch similarity approach against Brute Force search using the *Cemetery* dataset mentioned, and the results are shown in **Table 1**.

**Table 1: Sketch similarity evaluation by recall and runtime using Cemetery dataset across 3000 query polygons**

| Grid Dimension | Recall @10 | Recall @100 | Recall @1000 | Time (seconds) |
|----------------|------------|-------------|--------------|----------------|
| 25 × 25 (414)  | 0.09       | 0.15        | 0.28         | 54.7           |
| 30 × 30 (583)  | 0.17       | 0.28        | 0.44         | 123.2          |
| 50 × 50 (1536) | 0.27       | 0.37        | 0.52         | 189.7          |
| 100 × 100 (5878) | 0.48     | 0.58        | 0.70         | 1088.5         |

From the experiment, we found that the encoded vectors alone do not adequately capture the Jaccard similarity. This is reflected in the lower recall values. However, recall improves when the length of the encoded sketch vectors increases for higher grid dimensions, which leads to high runtime in similarity search compared to the Brute Force (BF) method. 

For example, it took about **1088 seconds** using the sketch produced by a 100×100 grid. However, the BF search took only **461 seconds**. Using a KD-tree with these vectors can make the runtime faster, but the recall was still low. Thus, we opted for a **filter-and-refine approach** for similarity-based query processing in two distinct phases with minhashing and LSH.

---

**Table 2: Impact of binary sketch on Recall and runtime across multiple datasets when Hash Length = 5**

| Dataset   | Grid Dimension | Recall @10 | Recall @100 | LSH Time (s) | BF Time (s) |
|-----------|----------------|------------|-------------|--------------|-------------|
| Cemetery  | 16 × 16        | 0.92       | 0.90        | 291.7        | 465.2       |
| Cemetery  | 30 × 30        | 0.90       | 0.86        | 206.3        | 465.2       |
| Sports    | 30 × 30        | 0.99       | 0.99        | 1772.9       | 1780.0      |
| Parks     | 30 × 30        | 0.98       | 0.96        | 4240.5       | 4876.5      |

---

## 1.1.2 LSH with Binary Sketch

Our findings indicated that using a **binary sketch** yielded unfavorable results in terms of runtime. As the polygons are centered, some similar polygons tend to share identical binary sketches. This is because partial overlaps with grid cells are marked as 1 even though the overlap may be 10% or 90%. 

Binary encoding is not capable of correctly encoding partial overlaps. This leads to the scenario where multiple similar polygons end up having the same locality-sensitive hashes. This increases the number of false positives in the query phase, thus increasing the number of candidates to query against, which in turn increases the runtime.  

The experimental results in **Table 2** show that reducing the size of the grid leads to an aggregation of polygons that share the same hash, consequently increasing both the recall and the runtime. In some cases, this exceeds the time required by the Brute Force approach.

