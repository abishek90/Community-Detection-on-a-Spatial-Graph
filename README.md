# Community-Detection-on-a-Spatial-Graph
This is a code to detect communities in a spatial random graph

There are two algorithms implemented - A GBG algorithm that is published in our paper (http://abishek90.github.io/CommDet.pdf), and an improved spectral synchronization 
algorithm that is not published, but empirically performs better.

To use the algorithm, run the `simulation_spectral.py' script. The parameters for this are a large integer N, a real number \lambda and 0 \leq b < a \leq 1. These parameters can be set in the main function at the end of the code script 'simulation_spectral.py'. These parameters indicates that there are roughly \lambda N nodes in the graph, half of which are in one community and the rest in another. The nodes have an uniform location label in the set [-\sqrt{N}/2,\sqrt{N}/2]^2. Conditional on the locations and communities, two nodes of the same community at a distance of \sqrt{\log(n} are connected with probability a and two nodes of opposite communities with distance between them smaller than \sqrt{\log(n} are connected with probability b. Nodes farther than \sqrt{\log(n}} are not connected. The Community Detection problem refers to how and when can one recover the communities, from a random sample of the graph and locations. 
