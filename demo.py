import networkx as nx
import numpy as np
from scipy.linalg import ldl, solve_triangular

def detect_symmetric_subgraphs(G, mu=-1.0, eps=1e-8, tau=1e-5, density_threshold=0.9):
    """
    DetectSymmetricSubgraphs(G, μ, ε=1e-8, τ=1e-5)
    Returns list of vertex lists (symmetric subgraphs).
    """
    A = nx.to_numpy_array(G, dtype=float)
    n = A.shape[0]
    
    # Shifted matrix
    B = A - mu * np.eye(n)
    
    # LDL^T factorization (scipy returns L unit lower, D block-diag, perm)
    L, D, perm = ldl(B, lower=True, hermitian=True)
    
    # Scan D for near-zero 1×1 / 2×2 blocks -> kernel seeds
    kernels = []
    i = 0
    while i < n:
        if i + 1 < n and abs(D[i, i+1]) > 1e-12:          # 2×2 block
            block = D[i:i+2, i:i+2]
            eigvals, eigvecs = np.linalg.eig(block)
            for j in range(2):
                if abs(eigvals[j]) < eps:
                    u = np.zeros(n)
                    u[i:i+2] = eigvecs[:, j]
                    kernels.append(u)
            i += 2
        else:                                              # 1×1 block
            if abs(D[i, i]) < eps:
                u = np.zeros(n)
                u[i] = 1.0
                kernels.append(u)
            i += 1
    
    # Recover approximate eigenvectors
    evecs = []
    for u_seed in kernels:
        try:
            v = solve_triangular(L.T, u_seed, lower=False)   # L^T v = u
            v = v / np.linalg.norm(v) if np.linalg.norm(v) > 1e-12 else v
        except:
            v = u_seed.copy()
        
        x = np.zeros(n)
        x[perm] = v  # undo permutation
        evecs.append(x)
    
    # Build auxiliary graph H from supports
    H = nx.Graph()
    for x in evecs:
        supp = [k for k in range(n) if abs(x[k]) > tau]
        if len(supp) >= 2:
            H.add_nodes_from(supp)
            for a in supp:
                for b in supp:
                    if a < b:
                        H.add_edge(a, b)
    
    subgraphs = [sorted(list(c)) for c in nx.connected_components(H) if len(c) >= 2]
    
    # Refinement (handles overlaps/duplicates via density)
    refined = []
    for S in subgraphs:
        induced = G.subgraph(S)
        dens = nx.density(induced)
        if dens >= density_threshold:                    # e.g. 0.9 for near-clique
            refined.append(S)
        else:
            # Future: max-clique split or community detection here
            refined.append(S)                            # placeholder for now
    return refined