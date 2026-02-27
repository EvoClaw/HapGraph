"""
HapGraph: Bayesian Population Admixture Graph Inference
using Joint F-statistics and IBD Sharing Statistics.

Estimates topology, admixture proportions (alpha), and admixture timing (T)
from phased VCF data.
"""

__version__ = "0.1.0-dev"
__author__  = "HapGraph Developers"

from .topology.nj_tree import nj_tree_from_f2
from .topology.greedy_search import greedy_admixture_search
from .inference.likelihood import HapGraphLikelihood
