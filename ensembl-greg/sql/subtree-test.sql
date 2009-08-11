SELECT
       ptnp.node_id AS p_node,
       leaf_count(ptnp.node_id) AS p_leaves,
       gappiness(ptnp.node_id) AS p_gaps,
       gappiness(ptna.node_id) AS a_gaps,
       gappiness(ptnb.node_id) AS b_gaps,
       avg_omega_for_node(ptnp.node_id,101) AS p_avg,
       avg_omega_for_node(ptna.node_id,101) AS a_avg,
       avg_omega_for_node(ptnb.node_id,101) AS b_avg,
       abs(avg_omega_for_node(ptna.node_id,101)-avg_omega_for_node(ptnb.node_id,101)) AS diff
  FROM
    node_set_member nsm, protein_tree_node ptnp, protein_tree_node ptna, protein_tree_node ptnb
  WHERE
    nsm.node_set_id=2 AND
    ptnp.node_id=nsm.node_id AND
    ptna.parent_id=ptnp.node_id AND
    ptnb.parent_id=ptnp.node_id AND
    ptna.node_id != ptnb.node_id AND
    ptna.node_id > ptnb.node_id AND
    avg_omega_for_node(ptnp.node_id,101) IS NOT NULL AND
    avg_omega_for_node(ptna.node_id,101) IS NOT NULL AND
    avg_omega_for_node(ptnb.node_id,101) IS NOT NULL
;
