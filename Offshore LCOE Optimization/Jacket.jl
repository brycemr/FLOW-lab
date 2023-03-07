function compute(num_legs)
    n_legs = num_legs
    material_name = "steel"
    eps = 1e-8

    # # Unpack inputs and ignore ghost nodes
    # c_trans = inputs["transition_piece_cost"]
    # leg_nodes = inputs["leg_nodes"]
    # connect_mat = np.int_(inputs["jacket_elem_N"][:-n_legs, :])
    # L = inputs["jacket_elem_L"][:-n_legs]
    # D = inputs["jacket_elem_D"][:-n_legs].copy()
    # t = inputs["jacket_elem_t"][:-n_legs]
    # xyz = inputs["jacket_nodes"]
    # m_total = inputs["jacket_mass"]
    # n_edges = L.size
    # N0 = connect_mat[:, 0]
    # N1 = connect_mat[:, 1]

    # # Compute costs based on "Optimum Design of Steel Structures" by Farkas and Jarmai
    # imat = discrete_inputs["material_names"].index(material_name)
    # k_m = inputs["unit_cost_mat"][imat]  # 1.1 # USD / kg carbon steel plate
    # k_f = inputs["labor_cost_rate"]  # 1.0 # USD / min labor
    # k_p = inputs["painting_cost_rate"]  # USD / m^2 painting
    # k_e = 0.064  # Industrial electricity rate $/kWh https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_5_6_a
    # e_f = 15.9  # Electricity usage kWh/kg for steel
    # # e_fo = 26.9  # Electricity usage kWh/kg for stainless steel
    # theta = 3  # Manufacturing difficulty factor

    # # Prep vectors for elements and legs for some vector math
    # edge_vec = xyz[N0, :] - xyz[N1, :]
    # leg_vec = np.squeeze(leg_nodes[:, -1, :] - leg_nodes[:, 0, :])
    # leg_L = np.linalg.norm(leg_vec, axis=1)

    # # Get the angle of intersection between all edges (as vectors) and all legs (as vectors)
    # vec_vals = np.dot(edge_vec, leg_vec.T) / np.outer(L, leg_L)
    # leg_alpha = np.arccos(np.minimum(np.maximum(vec_vals, -1.0), 1.0))
    # # If angle of intersection is close to 0 or 180, the edge is part of a leg
    # tol = np.deg2rad(5)
    # idx_leg = np.any((np.abs(leg_alpha) < tol) | (np.abs(leg_alpha - np.pi) < tol), axis=1)
    # D[idx_leg] = 0.0  # No double-counting time for leg elements since single piece

    # # Now determine which angle to use based on which leg a given edge node is on
    # edge_alpha = 0.5 * np.pi * np.ones(n_edges)
    # for k in range(n_legs):
    #     tol = 1e-2 * leg_L[k]
    #     sec1 = np.linalg.norm(np.squeeze(leg_nodes[k, -1, :])[np.newaxis, :] - xyz[N0, :], axis=1)
    #     sec2 = np.linalg.norm(xyz[N0, :] - np.squeeze(leg_nodes[k, 0, :])[np.newaxis, :], axis=1)
    #     on_leg_k = np.abs(leg_L[k] - sec1 - sec2) < tol
    #     edge_alpha[on_leg_k] = leg_alpha[on_leg_k, k]
    # edge_alpha = np.minimum(edge_alpha, np.pi - edge_alpha) + eps

    # # Run manufacturing time estimate functions
    # weld_L = 2 * np.pi * D / np.sin(edge_alpha)  # Multiply by 2 for both ends
    # n_pieces = n_edges - np.count_nonzero(D)
    # t_cut = 2 * manu.steel_tube_cutgrind_time(theta, 0.5 * D, t, edge_alpha)  # Multiply by 2 for both ends
    # t_weld = manu.steel_tube_welding_time(theta, n_pieces, m_total, weld_L, t)
    # t_manufacture = t_cut + t_weld
    # K_f = k_f * t_manufacture

    # # Cost step 5) Painting by surface area
    # theta_p = 2
    # K_p = k_p * theta_p * np.pi * np.sum(D[:-n_legs] * L[:-n_legs])

    # # Material cost with waste fraction, but without outfitting,
    # K_m = 1.21 * np.sum(k_m * m_total)

    # # Electricity usage
    # K_e = k_e * (e_f * m_total)  # + e_fo * (coeff - 1.0) * m_total

    # # Assemble all costs for now
    # tempSum = K_m + K_e + K_p + K_f

    # # Capital cost share from BLS MFP by NAICS
    # K_c = 0.118 * tempSum / (1.0 - 0.118)

    # outputs["jacket_cost"] = K_c + tempSum + c_trans
    # outputs["structural_cost"] = outputs["jacket_cost"] + inputs["tower_cost"]

end