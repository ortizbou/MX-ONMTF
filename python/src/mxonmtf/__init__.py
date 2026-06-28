"""MX-ONMTF: Community detection in multiplex networks.

Multiplex Orthogonal Nonnegative Matrix Tri-Factorization for identifying
common and private communities across multiplex network layers.

Reference:
    M. Ortiz-Bouza and S. Aviyente, "Community Detection in Multiplex Networks
    Based on Orthogonal Nonnegative Matrix Tri-Factorization," IEEE Access,
    vol. 12, pp. 6423-6436, 2024.
"""

from mxonmtf.metrics import modularity_density, normalized_mutual_information
from mxonmtf.utils import normalized_adjacency, row_normalize

__version__ = "0.1.0"
