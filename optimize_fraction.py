import numpy as np
from optimize import *
import gurobipy as gp
import math

### use gurobi

def create_sum_same_clone(tree_list, node_list, gene2idx, tree_freq_list=None):
    num_tree = len(tree_list)
    if tree_freq_list is None:
        tree_freq_list = np.ones(num_tree)
    sum_tree_freq = np.sum(tree_freq_list)
    for i in range(num_tree):
        tree = tree_list[i]
        node_dict = node_list[i]
        sam_clo_matrix = create_same_clone_matrix(tree, node_dict, gene2idx)
        #sam_clo_matrix /= (np.sum(sam_clo_matrix, axis=0) + 1) # normalization
        if i == 0:
            sam_clo_matrix_sum = np.zeros(sam_clo_matrix.shape)
        sam_clo_matrix_sum += sam_clo_matrix * tree_freq_list[i]
    return sam_clo_matrix_sum


def create_gene_fraction_array(tree, node_dict, clonal_freq, gene2idx, focus_sample_idx=0):
    gene_fraction_array = np.zeros((len(gene2idx.keys())))
    root = root_searching(tree)
    for node, muts_list in node_dict.items():
        if node != root:
            muts_idx = [gene2idx[mut] for mut in muts_list]
            print(clonal_freq[node][0][focus_sample_idx])
            gene_fraction_array[np.array(muts_idx)] = clonal_freq[node][0][focus_sample_idx]
    return gene_fraction_array


def create_concat_gene_fraction(tree_list, node_list, clonal_freq_list, gene2idx, tree_freq_list=None,focus_sample_idx=0):
    num_trees = len(tree_list)
    num_genes = len(gene2idx.keys())
    if tree_freq_list is None:
        tree_freq_list = np.ones(num_trees)
    sum_tree_freq = np.sum(tree_freq_list)
    gene_fraction_matrix = np.zeros((num_trees, num_genes))
    for i in range(num_trees):
        tree = tree_list[i]
        node_dict = node_list[i]
        clonal_freq = clonal_freq_list[i]
        gene_fraction_array = create_gene_fraction_array(tree, node_dict, clonal_freq, gene2idx, focus_sample_idx)
        gene_fraction_matrix[i] = gene_fraction_array
    return gene_fraction_matrix


def create_gene_variance_matrix(gene_fraction_matrix, read_depth=10000):
    return read_depth * gene_fraction_matrix * (1-gene_fraction_matrix)


def create_concat_relation_matrix(tree_list, node_list, gene2idx, tree_freq_list=None):
    num_tree = len(tree_list)
    num_gene = len(gene2idx.keys())
    if tree_freq_list is None:
        tree_freq_list = np.ones(num_tree)
    relation_matrix_full = np.zeros((num_tree, num_gene, num_gene))
    for i in range(num_tree):
        tree = tree_list[i]
        node_dict = node_list[i]
        relation_matrix_full[i, :, :] = create_ancestor_descendant_matrix(tree, node_dict, gene2idx)
    return relation_matrix_full

def optimize_tree_distribution(F, R,  n_genes, n_markers, read_depth, lam1, lam2, tree_freq_list, subset_list=None):
    model = gp.Model('opt_tree')
    V_sqr = create_gene_variance_matrix(F, read_depth)
    n_trees = F.shape[0]
    F_12 = F[:,:, np.newaxis]
    F_23 = np.transpose(F)[np.newaxis, :]
    #print(F_12.shape, F_23.shape)
    V_sqr_12 = V_sqr[:,:, np.newaxis]
    V_sqr_23 = np.transpose(V_sqr)[np.newaxis, :]
    frac_diff_matrix = F_12 - F_23
    #print(frac_diff_matrix.shape)
    var_sum_matrix = V_sqr_12 + V_sqr_23
    ratio_matrix = frac_diff_matrix ** 2 / var_sum_matrix
    log_likelihood_matrix = - 1 / 2 * read_depth ** 2 * ratio_matrix - 1 / 2 * np.log(var_sum_matrix)
    for i in range(n_trees):
        log_likelihood_matrix[i, :, i] = 0  # zero out the diagonal

    R_12 = R[np.newaxis, :, :, :]
    R_23 = R[:, np.newaxis, :, :]
    R_abs_diff = np.abs(R_12 - R_23)
    print(R_abs_diff.shape)
    z = get_gp_1d_arr_bin_var(model, n_genes)
    sum_struct = model.addVar(vtype=gp.GRB.INTEGER, lb=0, ub=n_trees**2*n_markers**2)
    Obj_frac = model.addVar(vtype=gp.GRB.CONTINUOUS)
    Obj_struct = model.addVar(vtype=gp.GRB.CONTINUOUS)
    model.addConstr(Obj_struct == sum_struct * np.log(10))
    model.addConstr(gp.quicksum([R_abs_diff[i, k, l, m]*z[l]*z[m]*tree_freq_list[i]*tree_freq_list[k] for l in range(n_genes)
                    for m in range(n_genes) for i in range(n_trees) for k in range(n_trees)])== sum_struct, name='obj_tree_struct_constraint')
    model.addConstr(- gp.quicksum([z[j] * log_likelihood_matrix[i, j, k] for i in range(n_trees)
                    for j in range(n_genes) for k in range(n_trees)]) == Obj_frac, name='obj_fraction_constraint')
    if subset_list is not None:
        subset_constraints(model, z, subset_list, n_markers, n_genes)
    else:
        model.addConstr(gp.quicksum([z[i] for i in range(n_genes)]) == n_markers, name='n_marker constraint')
    model.setObjective(gp.quicksum([lam1*Obj_frac, lam2*Obj_struct]), gp.GRB.MAXIMIZE)
    #model.setObjective(Obj_struct, gp.GRB.MAXIMIZE)
    model.optimize()
    return Obj_frac.X, Obj_struct.X, return_value_1d(z)


def select_markers_tree_gp(gene_list, n_markers, tree_list, node_list, clonal_freq_list, gene2idx, tree_freq_list,
                           read_depth=10000, lam1=0.001, lam2=1,focus_sample_idx=0, subset_list=None):
    F = create_concat_gene_fraction(tree_list, node_list, clonal_freq_list, gene2idx, tree_freq_list, focus_sample_idx)
    R = create_concat_relation_matrix(tree_list, node_list, gene2idx)
    n_genes = len(gene_list)
    best_obj_frac, best_obj_struct, best_z = optimize_tree_distribution(F, R, n_genes, n_markers, read_depth, lam1, lam2, tree_freq_list, subset_list)
    print(best_obj_frac, best_obj_struct, best_z)
    best_z = np.round(best_z).astype(int)
    selected_markers = []
    for idx in range(len(best_z)):
        if best_z[idx] == 1.0:
            selected_markers.append(gene_list[idx])
    return selected_markers, best_obj_frac, best_obj_struct


def optimize_fraction_weighted_single(E, M, F_hat, n_genes, n_markers):
    print(E, M, F_hat)
    k = E.shape[0]
    model = gp.Model('opt_frac')
    z = get_gp_1d_arr_bin_var(model, n_genes)
    x = get_gp_1d_arr_int_var(model, k)
    obj = model.addVar(vtype=gp.GRB.CONTINUOUS)
    x_bin = get_bin_1d(model, x, k)
    for i in range(k):
        model.addConstr(gp.quicksum([M[i, j] * z[j] for j in range(n_genes)]) == x[i])
    I = np.identity(k)
    I[0, 0] = 0      # set the normal clone to be 1 and no need to measure
    E_hat = E + I
    y = get_gp_1d_arr_int_var(model, k)
    t = get_gp_1d_arr_bin_var(model, k)
    for i in range(k):
        model.addConstr(gp.quicksum([E_hat[i, j] * (1 - x_bin[j]) for j in range(k)]) == y[i])
    y_bin = get_bin_1d(model, y, k)
    for i in range(k):
        model.addConstr(t[i] == 1 - y_bin[i])
    model.addConstr(gp.quicksum([z[i] for i in range(n_genes)]) == n_markers, name='n_marker constraint')
    model.addConstr(obj == gp.quicksum(t[i] * F_hat[i] for i in range(k)))
    model.setObjective(obj, gp.GRB.MAXIMIZE)
    model.optimize()
    return return_value_1d(z), obj.X


#TODO MODIFY THIS CONTENT OF THE FUNCTION
def calculate_fraction_weighted_single(z, E_list, M_list, F_hat_list, tree_freq_list, n_genes, n_markers, subset_list=None):
    n_trees = len(E_list)
    obj_list = np.zeros((n_trees))
    x_list = []
    x_bin_list = []
    y_list = []
    y_bin_list = []
    t_list = []
    for tree_idx, (E, M, F_hat, tree_freq) in enumerate(zip(E_list, M_list, F_hat_list, tree_freq_list)):
        k = E.shape[0]
        x = np.dot(M, z)
        x_bin = x.astype(bool)
        I = np.identity(k)
        I[0, 0] = 0  # set the normal clone to be 1 and no need to measure
        E_hat = E + I
        y = np.sum(E_hat * (1 - x_bin), axis=0)
        t = 1 - y.astype(bool)
        obj_list[tree_idx] = np.sum(np.dot(t * F_hat))
        x_list.append(x)
        x_bin_list.append(x_bin)
        y_list.append(y)
        t_list.append(t)
    final_mean_obj = np.dot(obj_list,tree_freq_list)
    return final_mean_obj/sum(tree_freq_list)
def calculate_fraction_weighted_overall(z, E_list, M_list, F_hat_list, tree_freq_list, n_genes, n_markers, subset_list=None):
    n_trees = len(E_list)
    obj_list = np.zeros((n_trees))
    x_list = []
    x_bin_list = []
    y_list = []
    y_bin_list = []
    t_list = []
    for tree_idx, (E, M, F_hat, tree_freq) in enumerate(zip(E_list, M_list, F_hat_list, tree_freq_list)):
        k = E.shape[0]
        x = np.dot(M, z)
        x_bin = x.astype(bool)
        I = np.identity(k)
        I[0, 0] = 0  # set the normal clone to be 1 and no need to measure
        E_hat = E + I
        y = np.sum(E_hat * (1 - x_bin), axis=0)
        t = 1 - y.astype(bool)
        obj_list[tree_idx] = np.sum(np.dot(t * F_hat))
        x_list.append(x)
        x_bin_list.append(x_bin)
        y_list.append(y)
        t_list.append(t)
    final_mean_obj = np.dot(obj_list,tree_freq_list)
    return final_mean_obj/sum(tree_freq_list)

def optimize_fraction_weighted_overall(E_list, M_list, F_hat_list, tree_freq_list, n_genes, n_markers, subset_list=None):
    model = gp.Model('opt_frac')
    z = get_gp_1d_arr_bin_var(model, n_genes)
    if subset_list is not None:
        subset_constraints(model, z, subset_list, n_markers, n_genes)
    else:
        model.addConstr(gp.quicksum([z[i] for i in range(n_genes)]) == n_markers, name='n_marker constraint')
    n_trees = len(E_list)
    obj_list = get_gp_1d_arr_cont_var(model, n_trees, 0)
    x_list = []
    x_bin_list = []
    y_list = []
    y_bin_list = []
    t_list = []
    for tree_idx, (E, M, F_hat, tree_freq) in enumerate(zip(E_list, M_list, F_hat_list, tree_freq_list)):
        k = E.shape[0]
        x = get_gp_1d_arr_int_var(model, k)
        x_bin = get_bin_1d(model, x, k)
        for i in range(k):
            model.addConstr(gp.quicksum([M[i, j] * z[j] for j in range(n_genes)]) == x[i])
        I = np.identity(k)
        I[0, 0] = 0  # set the normal clone to be 1 and no need to measure
        E_hat = E + I
        y = get_gp_1d_arr_int_var(model, k)
        t = get_gp_1d_arr_bin_var(model, k)
        for i in range(k):
            model.addConstr(gp.quicksum([E_hat[i, j] * (1 - x_bin[j]) for j in range(k)]) == y[i])
        y_bin = get_bin_1d(model, y, k)
        for i in range(k):
            model.addConstr(t[i] == 1 - y_bin[i])
        model.addConstr(gp.quicksum(t[i] * F_hat[i] for i in range(k)) == obj_list[tree_idx])
        x_list.append(x)
        x_bin_list.append(x_bin)
        y_list.append(y)
        y_bin_list.append(y_bin)
        t_list.append(t)
    # TO-DO: Update the objective function with weights based on how medically important the subclones are.
    model.setObjective(gp.quicksum(obj_list[i] * tree_freq_list[i] for i in range(n_trees)), gp.GRB.MAXIMIZE)
    model.update()
    model.optimize()
    return return_value_1d(z), return_value_1d(obj_list)

def calculate_markers_fractions_weighted_overall(gene_list, n_markers, tree_list, node_list, clonal_freq_list_new_gt, gene2idx, tree_freq_list, subset_list=None):
    k_list = create_k_list(node_list)
    E_list = tree2E_list(tree_list, k_list)
    n_genes = len(gene_list)
    M_list = create_M_list(node_list, gene2idx, n_genes)
    F_list, F_hat_list = create_F_F_hat_list_from_pcr(clonal_freq_list_new_gt, tree_list)
    final_mean_obj = calculate_fraction_weighted_overall(E_list, M_list, F_hat_list, tree_freq_list, n_genes, n_markers, subset_list)
    return final_mean_obj

def select_markers_fractions_weighted_overall(gene_list, n_markers, tree_list, node_list, clonal_freq_list, gene2idx, tree_freq_list, subset_list=None, sample_idx=0):
    k_list = create_k_list(node_list)
    E_list = tree2E_list(tree_list, k_list)
    n_genes = len(gene_list)
    M_list = create_M_list(node_list, gene2idx, n_genes)
    F_list, F_hat_list = create_F_F_hat_list(clonal_freq_list, tree_list, sample_idx)
    best_z, obj_list = optimize_fraction_weighted_overall(E_list, M_list, F_hat_list, tree_freq_list, n_genes, n_markers, subset_list)
    selected_markers = []
    for idx in range(len(best_z)):
        if best_z[idx] == 1:
            selected_markers.append(gene_list[idx])
    return selected_markers, np.mean(obj_list)

def select_markers_fractions_weighted_single(gene_list, n_markers, tree_list, node_list, clonal_freq_list, gene2idx, tree_freq_list, idx_best, sample_idx=0):
    k_list = create_k_list(node_list)
    E_list = tree2E_list(tree_list, k_list)
    n_genes = len(gene_list)
    M_list = create_M_list(node_list, gene2idx, n_genes)
    F, F_hat_list = create_F_F_hat_list(clonal_freq_list, tree_list, sample_idx)
    k, E, M, F_hat = k_list[idx_best], E_list[idx_best], M_list[idx_best], F_hat_list[idx_best]
    best_z, obj = optimize_fraction_weighted_single(E, M, F_hat, n_genes, n_markers)
    selected_markers = []
    for idx in range(len(best_z)):
        if best_z[idx] == 1:
            selected_markers.append(gene_list[idx])
    return selected_markers, obj


def subset_constraints(model, z, subset_list, n_markers, n_genes):
    model.addConstr(gp.quicksum([z[i] for i in subset_list]) == n_markers)
    for j in range(n_genes):
        if j not in subset_list:
            model.addConstr(z[j] == 0)


def optimize_fraction_unweighted_overall(S, n_genes, n_markers):
    model = gp.Model('opt_frac')
    z = get_gp_1d_arr_bin_var(model, n_genes)
    model.addConstr(gp.quicksum([z[i] for i in range(n_genes)]) == n_markers, name='n_marker constraint')
    model.setObjective(gp.quicksum(z[i] * S[i, j] * z[j] for i in range(n_genes) for j in range(n_genes)), gp.GRB.MINIMIZE)
    model.optimize()
    return return_value_1d(z)


def select_markers_fractions_gp(gene_list, n_markers, tree_list, node_list, gene2idx, tree_freq_list):
    S = create_sum_same_clone(tree_list, node_list, gene2idx, tree_freq_list)
    print(S)
    n_genes = len(gene_list)
    best_z = optimize_fraction_unweighted_overall(np.tril(S), n_genes, n_markers)
    #best_z = optimize_fraction(S, n_genes, n_markers)
    selected_markers = []
    for idx in range(len(best_z)):
        if best_z[idx] == 1:
            selected_markers.append(gene_list[idx])
    return selected_markers


###############  Utility functions  ##############
def get_gp_1d_arr_bin_var(model, m):
    X = np.empty((m), dtype=gp.Var)
    for i in range(m):
        X[i] = model.addVar(vtype=gp.GRB.BINARY)
    return X


def get_gp_1d_arr_cont_var(model, m, lower_bound,):
    X = np.empty((m), dtype=gp.Var)
    for i in range(m):
        X[i] = model.addVar(vtype=gp.GRB.CONTINUOUS, lb=lower_bound)
    return X


def get_gp_2d_arr_bin_var(model, m, n,):
    X = np.empty((m, n), dtype=gp.Var)
    for i in range(m):
        for j in range(n):
            X[i, j] = model.addVar(vtype=gp.GRB.BINARY)
    return X

def get_gp_3d_arr_bin_var(model, m, n, p,):
    X = np.empty((m, n, p), dtype=gp.Var)
    for i in range(m):
        for j in range(n):
            for k in range(p):
                X[i, j, k] = model.addVar(vtype=gp.GRB.BINARY)
    return X

def get_gp_4d_arr_bin_var(model, m, n, p, q):
    X = np.empty((m, n, p, q), dtype=gp.Var)
    for i in range(m):
        for j in range(n):
            for k in range(p):
                for l in range(q):
                    X[i, j, k, l] = model.addVar(vtype=gp.GRB.BINARY)
    return X


def get_gp_1d_arr_int_var(model, m,):
    X = np.empty((m), dtype=gp.Var)
    for i in range(m):
        X[i] = model.addVar(vtype=gp.GRB.INTEGER)
    return X

def get_gp_2d_arr_int_var(model, m, n,):
    X = np.empty((m, n), dtype=gp.Var)
    for i in range(m):
        for j in range(n):
            X[i, j] = model.addVar(vtype=gp.GRB.INTEGER)
    return X


def get_gp_3d_arr_int_var(model, m, n, p, vmin, vmax):
    X = np.empty((m, n, p), dtype=gp.Var)
    for i in range(m):
        for j in range(n):
            for k in range(p):
                if vmax is None:
                    X[i, j, k] = model.addVar(vtype=gp.GRB.INTEGER)
                else:
                    X[i, j, k] = model.addVar(lb=vmin, ub=vmax, vtype=gp.GRB.INTEGER)
    return X

def get_gp_4d_arr_int_var(model, m, n, p, q, vmin, vmax):
    X = np.empty((m, n, p, q), dtype=gp.Var)
    for i in range(m):
        for j in range(n):
            for k in range(p):
                for l in range(q):
                    if vmax is None:
                        X[i, j, k, l] = model.addVar(vtype=gp.GRB.INTEGER)
                    else:
                        X[i, j, k, l] = model.addVar(lb=vmin, ub=vmax, vtype=gp.GRB.INTEGER)
    return X

def get_abs(model, x,):
    x_abs = model.addVar(vtype=gp.GRB.CONTINUOUS)
    model.addConstr(x_abs, gp.GRB.GREATER_EQUAL, x)
    model.addConstr(x_abs, gp.GRB.GREATER_EQUAL, -1 * x)
    return x_abs

def get_bin_2d(model, X, vmax):
    m, n = X.shape
    Y = get_gp_2d_arr_bin_var(model, m, n)  # Y = 0 if X == 0. Y = 1 if X != 0
    num_bits = int(math.floor(math.log(vmax, 2))) + 1  # maximum number of bits required
    Z = get_gp_3d_arr_bin_var(model, m, n, num_bits)  # bit representation of X
    for i in range(m):
        for j in range(n):  # set Z as bit representation
            model.addConstr(gp.quicksum([Z[i, j, b] * 2 ** b for b in range(num_bits)]) == X[i, j])
            for b in range(num_bits):  # Y must be 1 if any bits are 1
                model.addConstr(Z[i, j, b] <= Y[i, j])  # Y must be 0 if all bits are 0
            model.addConstr(Y[i, j] <= gp.quicksum([Z[i, j, b] for b in range(num_bits)]))
    return Y

def get_bin_1d(model, X, vmax):
    m = X.shape[0]
    Y = get_gp_1d_arr_bin_var(model, m)  # Y = 0 if X == 0. Y = 1 if X != 0
    num_bits = int(math.floor(math.log(vmax, 2))) + 1  # maximum number of bits required
    Z = get_gp_2d_arr_bin_var(model, m, num_bits)  # bit representation of X
    for i in range(m):  # set Z as bit representation
        model.addConstr(gp.quicksum([Z[i, b] * 2 ** b for b in range(num_bits)]) == X[i])
        for b in range(num_bits):  # Y must be 1 if any bits are 1
            model.addConstr(Z[i, b] <= Y[i])  # Y must be 0 if all bits are 0
        model.addConstr(Y[i] <= gp.quicksum([Z[i, b] for b in range(num_bits)]))
    return Y


def return_value_1d(x):
    n = x.shape[0]
    z = np.empty(n)
    for i in range(n):
        z[i] = x[i].X
    return z

def return_value_2d(x):
    m, n = x.shape
    z = np.empty((m, n))
    for i in range(m):
        for j in range(n):
            z[i, j] = x[i, j].X
    return z

def return_value_3d(x):
    m, n, p = x.shape
    z = np.empty((m, n, p))
    for i in range(m):
        for j in range(n):
            for k in range(p):
                z[i, j, k] = x[i, j, k].X
    return z

def return_value_4d(x):
    m, n, p, q = x.shape
    z = np.empty((m, n, p, q))
    for i in range(m):
        for j in range(n):
            for k in range(p):
                for l in range(q):
                    z[i, j, k, l] = x[i, j, k, l].X
    return z

### try to use cvxpy to solve but failed
### keep getting errors:
### cvxpy.error.DCPError: Problem does not follow DCP rules. Specifically:
### The objective is not DCP, even though each sub-expression is.
### You are trying to minimize a function that is concave.

# import cvxpy as cp
#
#
# def solve_opt_frac(S, n_genes, n_markers):
#     ones = np.ones(n_genes)
#     A = np.random.random(S.shape)
#     z = cp.Variable(n_genes,)
#     prob = cp.Problem(cp.Minimize(cp.quad_form(z, S)),[cp.sum(z) == n_markers])
#     prob.solve()
#     return z.value
#
#
# def select_markers_fractions_cxv(gene_list, n_markers, tree_list, node_list, gene2idx, tree_freq_list):
#     S = create_sum_same_clone(tree_list, node_list, gene2idx, tree_freq_list)
#     rand_perturb = np.random.uniform(0, 10**(-5), size= S.shape)
#     S += rand_perturb + np.transpose(rand_perturb)
#     print(S)
#     n_genes = len(gene_list)
#     best_z = solve_opt_frac(S, n_genes, n_markers)
#     selected_markers = []
#     for idx in range(len(best_z)):
#         if best_z[idx] == 1:
#             selected_markers.append(gene_list[idx])
#     return selected_markers
#
#
