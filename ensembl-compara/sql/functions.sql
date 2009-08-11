DROP FUNCTION dubious_dup;
CREATE FUNCTION dubious_dup (n int(20))
RETURNS int(20)
RETURN (SELECT count(*) FROM protein_tree_tag where node_id=n AND tag="dubious_duplication" AND value="1");


DROP FUNCTION leaf_count;
CREATE FUNCTION leaf_count (n int(20))
RETURNS int(20)
RETURN (SELECT COUNT(*) 
   FROM protein_tree_node n1,protein_tree_node n2,protein_tree_member m 
   WHERE n2.left_index BETWEEN n1.left_index and n1.right_index
# GJ 2009-03-03 - until we redo the leftright-indexing, keep this in:
   AND n2.root_id=n1.root_id
   AND n1.node_id=n
   AND m.node_id=n2.node_id
);

