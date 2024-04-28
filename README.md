# Communities
Implementation of the paper

```
Yang, Jaewon, and Jure Leskovec. 
"Defining and evaluating network communities based on ground-truth." 
Knowledge and Information Systems 42.1 (2015): 181-213.
```

## Dataset
https://snap.stanford.edu/data/index.html
- DBLP
- LiveJournal

## Implemented quality functions

### Scoring functions based on internal connectivity
- Internal Density
- Edges inside
- Average Degree
- Fraction over median degree (FOMD)
- Triangle Participation Ratio (TPR)

### Scoring functions based on external connectivity
- Expansion
- Cut Ratio

### Scoring functions that combine internal and external connectivity
- Conductance
- Normalized Cut
- Maximum-ODF (Out Degree Fraction)
- Average-ODF
- Flake-ODF

### Scoring function based on a network model
- Modularity

### Community goodness metrics
- Separability
- Density
- Cohesiveness
- Clustering coefficient

## Discovering communities from a seed node
Implemented the algorithm to discover communities from a seed node based on the paper.

## References
- https://github.com/GiulioRossetti/partition_quality
- https://github.com/Lab41/Circulo
- https://github.com/ozanarkancan/community
