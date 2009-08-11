package Bio::EnsEMBL::Compara::RunnableDB::ProteinFeatureCollection;

use strict;
use Cwd;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

use Bio::AlignIO;

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

my $dba;
my $pta;
my $params;
my $node_id;

sub fetch_input {
    my $self = shift;

    # Load up the Compara DBAdaptor.
    $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
    $pta = $dba->get_ProteinTreeAdaptor;

    ### DEFAULT PARAMETERS ###
    $params->{'input_table'} = 'protein_tree_member';
    $params->{'output_table'} = 'aln_tag';
    $params->{'go_output_table'} = 'go_terms';
    $params->{'action'} = "features";  # Options: "features", "go", or any combination.
    #########################
    
    # Fetch parameters from the two possible locations. Input_id takes precedence!
    # (this utility method is from AlignmentProcess.pm)
    $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
    $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);

    #########################

    $self->check_if_exit_cleanly;

    #
    # Load the node_id.
    #
    $node_id = $params->{'protein_tree_id'};
    $node_id = $params->{'node_id'} if (!defined $node_id);

    #
    # Some last-minute adjustments based on retry counts or somesuch.
    #

    #
    # Think of reasons why we want to fail the job.
    #

}

my $jar_dir = "/nfs/users/nfs_g/gj1/src/uniprot";

sub run {
    my $self = shift;
    $self->check_if_exit_cleanly;
    $self->{'start_time'} = time()*1000;
    
    my $tmp = $self->worker_temp_directory;
    
    if ($params->{'action'} =~ "features")
    {
	my $cwd = cwd();
	chdir($tmp);
	# Prepare some config files for the Java execution.
	$self->_output_config_stuff;
	my $java = "/software/bin/java";
	system("${java} -jar ${jar_dir}/uniProtExtraction.jar ${node_id}");
	chdir($cwd);
    }
    
    if ($params->{'action'} =~ "go")
    {
	eval {
	    # Collect go terms for this protein tree.
	    my $go_terms = $self->collect_go_terms();
	}
    }

    $self->{'end_time'} = time()*1000;
}


sub _output_config_stuff {
    my $self = shift;

    my $input_table = $params->{'input_table'};
    my $output_table = $params->{'output_table'};

    open(OUT,">ensembl.properties");
    print OUT qq^
comparaConfig = mysql-compara2.cnf
anonConfig = mysql-ensembl.cnf
inputDb = gj1_compara_53
outputDb = gj1_compara_53
alignmentTable = $input_table
outputTable = $output_table
	^;
    close(OUT);

    open(OUT,">mysql-compara2.cnf");
    print OUT qq^
host = ens-research
user = ensadmin
password = ensembl
port = 3306
db = gj1_compara_53
	^;
    close(OUT);

    open(OUT,">mysql-ensembl.cnf");
    print OUT qq^
host = ensembldb.ensembl.org
user = anonymous
port = 5306
db = ensembl_compara_53
^;
    close(OUT);

    open(OUT,">sitewise_aln_tag.sql");
    my $str = qq^
	CREATE table if not exists aln_tag
	(
	 node_id int(10) unsigned NOT NULL,
	 aln_position int(10) unsigned NOT NULL,
	 member_id int(10) unsigned NOT NULL,
	 stable_id varchar(20) default NULL,
	 source varchar(20) default NULL,
	 tag varchar(16) NOT NULL default '',
	 value varchar(16) default NULL,
	 feature_position int(10) unsigned default NULL,
	 KEY node_id (node_id),
	 KEY tag (tag),
#	 FOREIGN KEY (node_id, aln_position) REFERENCES sitewise_aln(node_id, aln_position),
	 UNIQUE tag_aln_id (node_id,aln_position,tag)
	 );
    ^;
    print OUT $str;

    open(OUT,">uniprot.properties");
    print OUT qq^
#PROXY INFO
#Leave username password empty if no credientials required
#username=
#password=
#proxy.host=wwwcache.sanger.ac.uk
#proxy.port=3128
	^;
    close(OUT);
}

sub write_output {
    my $self = shift;

    my $tree = $pta->fetch_node_by_node_id($node_id);

    #
    # Store a tag for the runtime.
    #
    my $runtime = $self->{'end_time'} - $self->{'start_time'};
    print "RUNTIME: $runtime\n";
    $tree->store_tag("ProteinFeatureCollection_runtime",$runtime);
}



sub collect_go_terms
{
    my $self = shift;

    my $tree = $pta->fetch_node_by_node_id($node_id);
    my @leaves = @{$tree->get_all_leaves};

    # GJ 2009-02-09 : priority for grabbing GO annotations. Human > mouse > zebrafish > drosophila > c.elegans
#    my @species_leaves = grep {$_->taxon_id == 9606 || $_->taxon_id == 10090 ||
#			       $_->taxon_id == 7955 || $_->taxon_id == 7227 ||
#			   $_->taxon_id == 6239} @leaves;
#    @leaves = @species_leaves if (scalar(@species_leaves) >= 1);
    my @leaves = grep {$_->taxon_id==9606} @leaves;

    my %go_terms;
    my $taxon_id;

    foreach my $leaf (@leaves)
    {
	$taxon_id = $leaf->taxon_id;
	my $ts = $leaf->transcript;
	next if (!defined $ts);

	my $db_entries = $ts->get_all_DBLinks;

	# Grep out all the GO xref entries.
	my @keepers = grep {$_->dbname eq "GO"} @{$db_entries};
	foreach my $db_e (@keepers)
	{
	    $self->insert_go_term($node_id,$leaf,$db_e);
	    sleep(0.1);
	}
    }
    
    my @gos = keys(%go_terms);
    return \@gos;
}

sub insert_go_term
{
    my $self = shift;
    my ($node_id,$leaf,$db_e) = @_;

    my $cmd = "INSERT IGNORE INTO go_terms (node_id,member_id,stable_id,go_term) values (?,?,?,?);";
    print $cmd."\n";
    my $sth = $dba->dbc->prepare($cmd);
    $sth->execute($node_id,$leaf->dbID,$leaf->stable_id,$db_e->display_id);
}


1;
