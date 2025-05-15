# Dynamic Maintenance of Strong Connectivity under Bounded Edge Failures and Insertions

This repository contains the full C++ implementation of the algorithms described in this paper.

## 🖋 Thesis Abstract
In dynamic networks, failures such as link breakdowns are common. The work tackles the challenge of maintaining **strong connectivity information** in directed graphs under **edge failures and updates**, using efficient preprocessing and querying strategies. It builds on recent theoretical advances and includes an original extension for answering **pairwise strong connectivity queries** under both deletions and insertions.

---

## 📌 Key Features

- **Fault-Tolerant SCC Computation**  
  Preprocesses the input graph to allow efficient recomputation of SCCs after up to `k` edge deletions.

- **Dynamic Connectivity Queries**  
  Quickly answers whether any two nodes remain strongly connected even after simultaneous deletions and insertions.

- **Efficient Algorithms & Structures**  
  Implements techniques like:
  - `k`-Fault Tolerant Reachability Subgraphs (k-FTRS)
  - Heavy Path Decomposition
  - Path-restricted SCC computations

- **Optimized for Real-Time Scenarios**  
  Enables rapid responses in failure-prone networks such as distributed systems, transportation grids, and communication infrastructures.

---

## 🧠 Based On

- Kosaraju’s and Tarjan’s algorithms for SCCs (used as benchmarks)
- k-FTRS framework by Baswana et al.
- Recent theoretical results on SCCs under failures with preprocessing time `O(2^k n log² n)` and data structure size `O(2^k n²)`
- My thesis contribution: **Pairwise strong connectivity query handling under both insertions and deletions**

---


