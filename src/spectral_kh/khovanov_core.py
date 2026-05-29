import numpy as np
import itertools
import math
import time
from collections import defaultdict
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import splu

# =================================================================
# CLASS: Smoothing
# Handles circle finding and basis generation for a specific state
# =================================================================
class Smoothing:
    """
    Represents a single resolution of the knot diagram.
    
    Attributes:
        state (tuple): The bit-string of 0s and 1s.
        circles (list): List of dictionaries containing circle 'min' and 'arcs'.
    """
    
    def __init__(self, state, pd_data):
        # Convert state to a standard tuple of integers
        if isinstance(state, str):
            self.state = tuple(int(bit) for bit in state)
        else:
            self.state = tuple(state)
            
        self.r = sum(self.state)  # Homological degree
        self.circles = self._compute_circles(pd_data)
        self.num_circles = len(self.circles)
        # Useful for matching c[m] in Mathematica
        self.min_labels = sorted([c['min'] for c in self.circles])

    def _compute_circles(self, pd_data):
        arcs = set()
        for x in pd_data: arcs.update(x)
        
        parent = {arc: arc for arc in arcs}
        def find(i):
            if parent[i] == i: return i
            parent[i] = find(parent[i]) # Path compression for efficiency
            return parent[i]

        def union(i, j):
            root_i, root_j = find(i), find(j)
            if root_i != root_j: parent[root_i] = root_j

        # Apply smoothings based on state bits
        for bit, (i, j, k, l) in zip(self.state, pd_data):
            if bit == 0: # 0-smoothing connects (i,j) and (k,l)
                union(i, j); union(k, l)
            else:        # 1-smoothing connects (i,l) and (j,k)
                union(i, l); union(j, k)

        circles_map = {}
        for arc in arcs:
            root = find(arc)
            circles_map.setdefault(root, []).append(arc)

        return sorted([{'min': min(c), 'arcs': sorted(c)} for c in circles_map.values()], key=lambda x: x['min'])

    def _add_edge(self, adj, u, v):
        """Helper to create bidirectional connections."""
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)

    def __repr__(self):
        return f"Smoothing(state={self.state}, circles={self.num_circles})"

    def get_circle_for_arc(self, arc_index):
        """
        Given an arc index (from the PD), finds the 'min' label 
        of the circle that contains it.
        """
        for circle in self.circles:
            if arc_index in circle['arcs']:
                return circle['min']
        return None  # Should not happen if PD is valid
    
    def generate_basis(self):
        """
        Generates the vector space basis using {1, X}.
        Returns a list of dictionaries.
        """
        basis = []
        # We now use '1' and 'X' as our labels
        options = [1, 'X']
        
        for labels in itertools.product(options, repeat=self.num_circles):
            # Calculate weight: 1 counts as +1, X counts as -1
            weight = sum(1 if L == 1 else -1 for L in labels)
            
            # Quantum degree: weight + homological degree (r)
            q_val = weight + self.r
            
            basis.append({
                'labels': labels, 
                'q': q_val,
                'state': self.state
            })
        return basis


# =================================================================
# CLASS: KhovanovEngine
# The main driver for loading PD data and computing the complex
# =================================================================
class KhovanovEngine:
    """
    The main driver for knot PD analysis and Khovanov Chain Complex generation.
    
    This class manages the 'Cube of Resolutions', calculates the grading shifts 
    based on crossing orientations, and provides the differential maps between 
    basis elements.
    """
    
    # -----------------------------------------------------------------
    # Initialization
    # -----------------------------------------------------------------
    def __init__(self, pd_list):
        self.pd = pd_list
        self.n = len(pd_list)
        
        self.crossing_signs = self._compute_crossing_signs()
        self.n_pos = self._count_pos()
        self.n_neg = self._count_neg()
        
        # Build the complex directly into chain_groups
        self.build_complex_by_degree()

    # --- Initialization Helpers ---
    
    def _compute_crossing_signs(self):
        """Return crossing signs using Spherogram's PD orientation when available."""
        try:
            from spherogram import Link

            link = Link(self.pd, check_planarity=False)
            signs = [int(c.sign) for c in link.crossings]
            if len(signs) == self.n and all(sign in (-1, 1) for sign in signs):
                return signs
        except Exception:
            pass

        signs = []
        for i, j, k, l in self.pd:
            if (j - l == 1) or (l - j > 1):
                signs.append(1)
            elif (l - j == 1) or (j - l > 1):
                signs.append(-1)
            else:
                raise ValueError(f"Could not determine crossing sign for {(i, j, k, l)}")
        return signs

    def _count_pos(self):
        """Counts positive crossings."""
        return sum(1 for sign in self.crossing_signs if sign == 1)

    def _count_neg(self):
        """Counts negative crossings."""
        return sum(1 for sign in self.crossing_signs if sign == -1)

    # -----------------------------------------------------------------
    # Smoothing Logic
    # -----------------------------------------------------------------
    def get_smoothing(self, state):
        """Creates a Smoothing object for a specific bit-string.""" 
        return Smoothing(state, self.pd)

    def _build_chain_complex(self):
        """Groups all 2^n resolutions by their raw homological degree r."""
        complex_dict = {r: [] for r in range(self.n + 1)}
        
        # itertools.product([0, 1], repeat=n) gives all vertices of the cube
        for state in itertools.product([0, 1], repeat=self.n):
            s_obj = Smoothing(state, self.pd)    # Creates a smoothing instance for each resolution
            complex_dict[s_obj.r].append(s_obj)  # Stores each smoothing in the corresponding homological degree r
        return complex_dict

    def build_complex_by_degree(self):
        """
        Populates a dictionary of all basis elements grouped by (h, q).
        This must be called before accessing self.chain_groups.
        """
        self.chain_groups = {}

        # 1. Iterate through every possible state in the Cube (2^n)
        for state in itertools.product([0, 1], repeat=self.n):
            s_obj = self.get_smoothing(state)
            basis_elements = s_obj.generate_basis()
            
            for b in basis_elements:
                # 2. Apply the Grading Shifts
                # h = r - n_neg
                # q = raw_q + n_pos - 2*n_neg
                h_shifted = s_obj.r - self.n_neg
                q_shifted = b['q'] + self.n_pos - 2 * self.n_neg
                
                degree_key = (h_shifted, q_shifted)
                
                # 3. Initialize the bucket if it doesn't exist
                if degree_key not in self.chain_groups:
                    self.chain_groups[degree_key] = []
                
                # Attach the circle names for the differential logic later
                b['min_labels'] = s_obj.min_labels
                self.chain_groups[degree_key].append(b)
        
        return self.chain_groups
        
    def get_active_bidegrees(self):
        """Returns the list of (h, q) tuples that have non-empty chain groups."""
        return sorted(self.chain_groups.keys())
    # -----------------------------------------------------------------
    # Differential Logic
    # -----------------------------------------------------------------    
    def get_differential(self, element: dict):
        """
        Calculates the complete differential d(element).
        Iterates through all possible 0->1 flips and applies signs.
        
        Returns:
            list[dict]: A list of resulting basis elements with 'coeff' keys.
        """
        results = []
        state = element['state']
        
        for k in range(self.n):
            if state[k] == 0:
                # 1. Calculate the Sign: (-1)^(number of 1s before index k)
                sign = (-1) ** sum(state[:k])
                
                # 2. Get algebra map results (Merge or Split)
                contributions = self._get_differential_contribution(element, k)
                
                # 3. Apply the sign to the results
                for term in contributions:
                    term['coeff'] = sign
                    results.append(term)
        return results
    
    def _get_differential_contribution(self, element: dict, bit_to_flip: int):
        """
        Determines if a bit-flip is a Merge or a Split and applies logic.
        
        Args:
            element (dict): A basis element dictionary containing:
                - 'state': tuple of 0s and 1s (the source vertex in the cube)
                - 'labels': tuple of 1s and 'X's (the algebra state of the circles)
                - 'min_labels': list of circle names in this smoothing
            bit_to_flip (int): The index of the crossing in self.pd to change 
                               from a 0-smoothing to a 1-smoothing.
                               
        Returns:
            list[dict]: A list of resulting basis elements in the next homological degree.
                        (Empty if the result is 0, or multiple elements if a Split occurs).
        """
        # Get the current smoothing and the neighbor smoothing
        s0 = self.get_smoothing(element['state'])
        new_state = list(element['state'])
        new_state[bit_to_flip] = 1    # flips the bit specified by bit_to_flip from 0 to 1
        s1 = self.get_smoothing(tuple(new_state))
        
        # Identify the arcs involved in the flip
        # Crossing X(i,j,k,l) at bit_to_flip
        i, j, k, l = self.pd[bit_to_flip]
        
        # In state 0, arc i and l are NOT connected. In state 1, they ARE.
        # Find which circles in s0 contain arcs i and l
        c_i = s0.get_circle_for_arc(i)
        c_l = s0.get_circle_for_arc(l)
        
        if c_i == c_l:
            # It's a SPLIT (One circle in s0 becomes two in s1)
            return self._apply_split(element, s0, s1, c_i)
        else:
            # It's a MERGE (Two circles in s0 become one in s1)
            return self._apply_merge(element, s0, s1, c_i, c_l)


    #---- Define Frobenius m map  -----
    def _apply_merge(self, element: dict, s0, s1, id_a, id_b):
        """
        Implements m: V ⊗ V -> V
        Rules: 1⊗1=1, 1⊗X=X, X⊗1=X, X⊗X=0
        """
        # 1. Identify which position in the labels tuple corresponds to our circles
        idx_a = s0.min_labels.index(id_a)
        idx_b = s0.min_labels.index(id_b)
        
        label_a = element['labels'][idx_a]
        label_b = element['labels'][idx_b]
        
        # 2. Apply the multiplication rule
        if label_a == 1 and label_b == 1:
            res_label = 1
        elif (label_a == 1 and label_b == 'X') or (label_a == 'X' and label_b == 1):
            res_label = 'X'
        else:
            return [] # X ⊗ X = 0 (The map vanishes)

        # 3. Build the new labels tuple for the target smoothing
        new_labels = []
        # We keep all labels from circles NOT involved in the merge
        # and add our new_label for the merged circle.
        for m_label in s1.min_labels:
            if m_label == min(id_a, id_b): # The "name" of the merged circle
                new_labels.append(res_label)
            else:
                # Find the original label for this untouched circle
                orig_idx = s0.min_labels.index(m_label)
                new_labels.append(element['labels'][orig_idx])
        
        return [{'labels': tuple(new_labels), 'state': s1.state}]

    
    #---- Define Frobenius comultiplication map  -----
    def _apply_split(self, element: dict, s0, s1, id_parent):
        """
        Implements Delta: V -> V ⊗ V
        Rules: 1 -> 1⊗X + X⊗1, X -> X⊗X
        """
        idx_parent = s0.min_labels.index(id_parent)
        label_parent = element['labels'][idx_parent]
        
        # Split results in one or two new basis elements
        results = []
        
        # Find the IDs of the two new circles in s1 that came from the parent
        # These are the IDs in s1.min_labels that weren't in s0.min_labels (plus the ID that stayed)
        child_ids = [m for m in s1.min_labels if m not in s0.min_labels or m == id_parent]
        # (Note: In a true split, s1 will have one more min_label than s0)

        # Define the outcome labels based on 1 -> 1⊗X + X⊗1 or X -> X⊗X
        outcomes = []
        if label_parent == 1:
            outcomes = [(1, 'X'), ('X', 1)]
        else:
            outcomes = [('X', 'X')]

        for left, right in outcomes:
            new_labels = []
            for m_label in s1.min_labels:
                if m_label == child_ids[0]:
                    new_labels.append(left)
                elif m_label == child_ids[1]:
                    new_labels.append(right)
                else:
                    orig_idx = s0.min_labels.index(m_label)
                    new_labels.append(element['labels'][orig_idx])
            results.append({'labels': tuple(new_labels), 'state': s1.state})
            
        return results

    # --- Matrix Constructor ---
    def get_differential_matrix(self, h, q):
        """
        Constructs a SPARSE matrix representation of the differential.
        Much faster for knots with 14+ crossings.
        """
        from scipy import sparse

        source_basis = self.chain_groups.get((h, q), [])
        target_basis = self.chain_groups.get((h + 1, q), [])

        if not source_basis or not target_basis:
            return sparse.csc_matrix((len(target_basis), len(source_basis)), dtype=int)

        target_lookup = {(tuple(b['state']), tuple(b['labels'])): i for i, b in enumerate(target_basis)}

        # Use Coordinate (COO) format to build the matrix efficiently
        rows = []
        cols = []
        data = []

        for col_idx, element in enumerate(source_basis):
            results = self.get_differential(element)
            for res in results:
                row_key = (tuple(res['state']), tuple(res['labels']))
                if row_key in target_lookup:
                    row_idx = target_lookup[row_key]
                    rows.append(row_idx)
                    cols.append(col_idx)
                    data.append(res['coeff'])

        # Convert to CSC for fast arithmetic
        matrix = sparse.coo_matrix((data, (rows, cols)), 
                                  shape=(len(target_basis), len(source_basis)), 
                                  dtype=int).tocsc()
        return matrix
    # -----------------------------------------------------------------
    # Compute Khovanov Poincare Polynomial and Jones polynomial
    # -----------------------------------------------------------------  
    def compute_khovanov_polynomial(self):
        """
        Calculates homology using the Laplacian spectrum (harmonic representatives)
        rather than rank-nullity on differentials for improved numerical stability.
        """
        betti_numbers = {}
        
        # Get all unique h and q values present in the complex
        all_h = sorted(list(set(h for h, q in self.chain_groups.keys())))
        all_q = sorted(list(set(q for h, q in self.chain_groups.keys())))
        
        for q in all_q:
            for h in all_h:
                # Skip if the chain group is empty
                dim = self.get_dim(h, q)
                if dim == 0:
                    continue
                    
                # Use your adaptive spectrum logic (Sonar Pinger)
                # 1e-7 is a very safe threshold given your 10^-12 gap
                betti, _ = self.get_low_spectrum_adaptive(h, q, tol=1e-7)
                
                if betti > 0:
                    betti_numbers[(h, q)] = betti
        
        # Use the helper function to turn the dict into the string
        return format_kh_poly(betti_numbers)

        return self._format_polynomial(betti_numbers)

    def _format_polynomial(self, betti_dict):
        """Formats the (h, q) results into a t, q polynomial string."""
        terms = []
        # Sort by h then q for a clean output
        for (h, q) in sorted(betti_dict.keys()):
            coeff = betti_dict[(h, q)]
            coeff_str = f"{coeff}" if coeff > 1 else ""
            
            t_str = f"t^{{{h}}}" if h != 0 else ""
            q_str = f"q^{{{q}}}" if q != 0 else ""
            
            # Combine terms logically
            term = f"{coeff_str}{t_str}{q_str}"
            if term == "": term = f"{coeff}" # Case for t^0 q^0
            terms.append(term)
            
        return " + ".join(terms)

    # -- Compute Jones polynomial ----
    def compute_jones_polynomial(self):
        """
        Collapses the Khovanov Polynomial into the Jones Polynomial 
        by setting t = -1.
        """
        jones_coeffs = {}
        
        # 1. Identify all active gradings
        all_h = sorted(list(set(h for h, q in self.chain_groups.keys())))
        all_q = sorted(list(set(q for h, q in self.chain_groups.keys())))

        for q in all_q:
            q_coeff = 0
            for h in all_h:
                # 2. Retrieve matrices inside the loop where h and q are defined
                m_curr = self.get_differential_matrix(h, q)
                m_prev = self.get_differential_matrix(h - 1, q)
                
                # Check matrix existence and identify dimension
                if m_curr is None:
                    continue
                
                dim_v = m_curr.shape[1]
                
                # 3. Perform rank calculations with correct 12-space indentation
                # Use .toarray() for sparse matrices to ensure matrix_rank works
                rank_curr = np.linalg.matrix_rank(m_curr.toarray() if sparse.issparse(m_curr) else m_curr, tol=1e-7) if m_curr.size > 0 else 0
                rank_prev = np.linalg.matrix_rank(m_prev.toarray() if sparse.issparse(m_prev) else m_prev, tol=1e-7) if m_prev.size > 0 else 0
                
                # 4. Calculate Homology Dimension: Dim(H) = (Dim(V) - Rank(d_curr)) - Rank(d_prev)
                kh_dim = (dim_v - rank_curr) - rank_prev
                
                # 5. Apply the t = -1 rule: (-1)^h * dim
                q_coeff += ((-1)**h) * kh_dim
            
            if q_coeff != 0:
                jones_coeffs[q] = q_coeff

        return self._format_jones(jones_coeffs)

    # -- Format Jones polynomial ----
    def _format_jones(self, jones_dict):
        """Formats the q-results into a standard polynomial string."""
        terms = []
        for q in sorted(jones_dict.keys()):
            coeff = jones_dict[q]
            sign = " + " if coeff > 0 else " - "
            abs_coeff = abs(coeff)
            
            coeff_str = f"{abs_coeff}" if abs_coeff != 1 or q == 0 else ""
            q_str = f"q^{{{q}}}" if q != 0 else ""
            
            terms.append(f"{sign}{coeff_str}{q_str}")
        
        return "".join(terms).strip(" + ")

    # -----------------------------------------------------------------
    # Compute Laplacians
    # -----------------------------------------------------------------  
    def get_upward_laplacian(self, h, q):
        m_out = self.get_differential_matrix(h, q)
        if m_out.shape[1] == 0:
            return sparse.csc_matrix((0, 0))
        return m_out.T @ m_out

    def get_downward_laplacian(self, h, q):
        m_in = self.get_differential_matrix(h - 1, q)
        if m_in.shape[0] == 0:
            # We need a zero matrix the size of the CURRENT space
            dim = len(self.chain_groups.get((h, q), []))
            return sparse.csc_matrix((dim, dim))
        return m_in @ m_in.T

    def compute_laplacian(self, h, q):
        # Keeps them sparse for efficient addition
        L_up = self.get_upward_laplacian(h, q)
        L_down = self.get_downward_laplacian(h, q)
        return L_up + L_down
    
    # ---------------------------------------------------------
    # -- Calculate spectrum of Laplacian -- 
    # OLD Version that truncated to 60 eigenvalues
    #def get_spectrum(self, h, q):
    #    laplacian = self.compute_laplacian(h, q)
    #    
    #    # Ensure we are treating it as a sparse matrix if it's large
    #    if laplacian.shape[0] > 100:
    #        try:
    #            # 'SM' finds the smallest magnitude eigenvalues (near zero)
    #            # which are the most important for homology
    #            evs = eigsh(laplacian, k=min(laplacian.shape[0]-1, 50), 
    #                        which='SM', return_eigenvectors=False)
    #        except:
    #            evs = np.linalg.eigvalsh(laplacian.toarray() if hasattr(laplacian, "toarray") else laplacian)
    #    else:
    #        evs = np.linalg.eigvalsh(laplacian.toarray() if hasattr(laplacian, "toarray") else laplacian)
    #            
    #    return evs
        
    def get_spectrum(self, h, q):
        """
        Computes the full spectrum of the Laplacian for bidegree (h, q).
        Uses a dense solver to ensure no eigenvalues are truncated, 
        which is critical for Reidemeister torsion calculations.
        """
        laplacian = self.compute_laplacian(h, q)
        dim = laplacian.shape[0]
        
        if dim == 0:
            return np.array([], dtype=float)
    
        # Convert sparse Laplacian to dense array and compute ALL eigenvalues.
        # np.linalg.eigvalsh is optimized for symmetric/Hermitian matrices.
        # This requires O(dim^2) memory and O(dim^3) time.
        return np.linalg.eigvalsh(laplacian.toarray())
    
    def get_low_spectrum_adaptive(self, h, q, tol=1e-9):
        laplacian = self.compute_laplacian(h, q)
        dim = laplacian.shape[0]
        
        if dim == 0:
            return 0, []
            
        # 1. Handle small matrices directly with a dense solver
        # Dense solvers are extremely fast for dim < 100 and avoid eigsh edge cases
        if dim < 100:
            evs = np.linalg.eigvalsh(laplacian.toarray())
            evs = np.sort(evs)
            betti = np.sum(np.abs(evs) < tol)
            gaps = evs[evs >= tol]
            return betti, gaps[:2] if len(gaps) > 0 else []

        # 2. Adaptive sparse solver for large matrices
        # 2. Adaptive sparse solver for large matrices
        k_search = 10
        while k_search < dim:
            # Use shift-invert (sigma) instead of which='SM' for massive speedup on Laplacians
            evs = eigsh(laplacian, k=k_search, sigma=-1e-5, return_eigenvectors=False)
            evs = np.sort(evs)
            
            betti = np.sum(np.abs(evs) < tol)
            
            # Check if we have 'cleared' the zero-kernel
            if np.abs(evs[-1]) > tol:
                next_two = evs[evs >= tol][:2]
                return betti, next_two
            
            # Prevent infinite loop if the Betti number is very large
            if k_search == dim - 1:
                break
                
            k_search = min(dim - 1, k_search * 2)
            
        # 3. Fallback if the sparse solver loop maxes out
        evs = np.linalg.eigvalsh(laplacian.toarray())
        evs = np.sort(evs)
        betti = np.sum(np.abs(evs) < tol)
        gaps = evs[evs >= tol]
        return betti, gaps[:2] if len(gaps) > 0 else []
            
    # -----------------------------------------------------------------
    # Torsion & Homology Analysis
    # -----------------------------------------------------------------
    def compute_rtorsion_report(self, tolerance=1e-7):
        """
        Calculates the Complex R-torsion tau^C(C) for each quantum grading q.
        Following Eq (5.1): Product of det'(Delta_i)^((-1)^{i+1} * i/2)
        """
        report = {}
        all_q = sorted(list(set(q for h, q in self.chain_groups.keys())))
        
        for q in all_q:
            total_log_val = 0.0
            total_betti = 0
            active_h = [h for h, q_val in self.chain_groups.keys() if q_val == q]
            
            for h in active_h:
                spec = self.get_spectrum(h, q)
                nonzero = spec[spec > tolerance]
                total_betti += len(spec) - len(nonzero)
                
                if len(nonzero) > 0:
                    # IMPLEMENTATION OF EQ (5.1) from Harmonic Khovanov Notes
                    exponent = ((-1)**(h + 1)) * (h / 2.0)
                    total_log_val += exponent * np.sum(np.log(nonzero))
            
            report[q] = {
                'tau_c': np.exp(total_log_val),
                'betti': total_betti,
                'is_acyclic': (total_betti == 0)
            }
        return report

    def get_homology_ranks(self, tolerance=1e-7):
        """
        Helper for students to see the Betti numbers in a grid format.
        """
        ranks = {}
        for (h, q) in self.get_active_bidegrees():
            spec = self.get_spectrum(h, q)
            rank = sum(1 for ev in spec if ev <= tolerance)
            if rank > 0:
                ranks[(h, q)] = rank
        return ranks
    # -----------------------------------------------------------------
    # Diagram invariants for comparisong to spectral gap
    # -----------------------------------------------------------------
    def get_turaev_genus(self):
        """
        Calculates the Turaev genus: g_T = (n + 2 - c0 - c1) / 2
        where c0 and c1 are the number of circles in the all-0 and all-1 resolutions.
        """
        # Create the all-0 and all-1 states as bit-tuples
        state_0 = tuple([0] * self.n)
        state_1 = tuple([1] * self.n)
        
        # Build the smoothings for those specific states and count their circles
        c0 = self.get_smoothing(state_0).num_circles
        cn = self.get_smoothing(state_1).num_circles
        
        turaev_genus = (self.n + 2 - c0 - cn) / 2
        return turaev_genus
    # -----------------------------------------------------------------
    # Info
    # -----------------------------------------------------------------   
    def get_info(self):
        """
        Returns a summary of the knot's properties and the grading 
        parameters used for the Khovanov complex.
        """
        # The total quantum shift used in: q_true = q_raw + r + n_pos - 2*n_neg
        quantum_shift = self.n_pos - 2 * self.n_neg
        
        return {
            "Total Crossings (n)": self.n,
            "Positive Crossings (n+)": self.n_pos,
            "Negative Crossings (n-)": self.n_neg,
            "Homological Shift (-n-)": -self.n_neg,
            "Quantum Shift (n+ - 2n-)": quantum_shift,
            "Homological Range": (0 - self.n_neg, self.n - self.n_neg),
            # The state space size (number of resolutions)
            "Complex Size": 2**self.n 
        }



# 1. The Wrapper Function (Calculates math on a single core)
def compute_bidegree_spectrum(args):
    engine, h, q = args
    
    # Internal timing for this specific bidegree
    start = time.perf_counter()
    spectrum = engine.get_spectrum(h, q)
    duration = time.perf_counter() - start
    
    return (h, q, sorted(spectrum, reverse=True), duration)

if __name__ == '__main__':
    # --- SETUP PHASE ---
    setup_start = time.perf_counter()
    
    figure_eight_pd = [(4,2,5,1), (8,6,1,5), (6,3,7,4), (2,7,3,8)]
    engine = KhovanovEngine(figure_eight_pd)
    
    all_keys = sorted(engine.chain_groups.keys(), key=lambda x: (x[1], x[0]))
    tasks = [(engine, h, q) for h, q in all_keys]
    
    setup_duration = time.perf_counter() - setup_start

    # --- PARALLEL SOLVE PHASE ---
    solve_start = time.perf_counter()
    
    with Pool() as pool:
        results = pool.map(compute_bidegree_spectrum, tasks)
        
    solve_duration = time.perf_counter() - solve_start

    # --- OUTPUT ---
    print(f"--- Performance Report ---")
    print(f"Setup Time (Basis Build): {setup_duration:.4f}s")
    print(f"Parallel Solve Time:      {solve_duration:.4f}s")
    print(f"Total Time:               {setup_duration + solve_duration:.4f}s\n")

    print(f"{'h':>3} | {'q':>3} | {'Time (s)':>10} | {'Spectrum (Decreasing)'}")
    print("-" * 65)
    
    for h, q, decreasing_spectrum, t_indiv in results:
        formatted_vals = [f"{val:.4f}" for val in decreasing_spectrum]
        list_str = "{" + ", ".join(formatted_vals) + "}"
        print(f"{h:3d} | {q:3d} | {t_indiv:10.4f} | {list_str}")



# =================================================================
# MULTIPROCESSING WORKERS (Used by driver.py)
# =================================================================

_worker_engine = None

def initialize_worker(pd_code):
    global _worker_engine
    _worker_engine = KhovanovEngine(pd_code)

def compute_single_spectrum_worker(bidegree):
    h, q = bidegree
    # FIX: Using the correct method name 'compute_laplacian'
    L = _worker_engine.compute_laplacian(h, q)
    
    # L is a sparse matrix here
    dim = L.shape[0]
    if dim == 0:
        return h, q, np.array([])

    try:
        # FOR REIDEMEISTER TORSION:
        # We need the full spectrum. Dense is safest for 14-crossings.
        if dim < 6000:
            # Densify only at the solve step to save memory during construction
            eigenvalues = np.linalg.eigvalsh(L.toarray())
        else:
            # Fallback for massive bidegrees
            from scipy.sparse.linalg import eigsh
            eigenvalues = eigsh(L, k=dim-1, which='SM', tol=1e-12)[0]
    except Exception as e:
        return h, q, np.array([-1.0]) 

    return h, q, eigenvalues
