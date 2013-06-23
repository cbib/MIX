#!/usr/bin/perl -w                                      -*- Perl -*-
#
# $Id$
#
# 

use Getopt::Long;
use List::Util qw(max min);
use List::MoreUtils qw(uniq);

my $do_help = 0;

my $blast_file = "";
my $pirsf_file = "";

GetOptions('help!'    => \$do_help,
           'blast=s'  => \$blast_file,
    );
 
if ( $do_help or $blast_file eq "") {
    print "usage: ",$0," --blast FILE [--help]\n";
    exit;
}

# Read the alignements and store them
# The BLAST file has the format 
# (Query id, Subject id, %identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score)
print STDERR "Reading BLAST alignements from $blast_file...";
open BLAST, "<$blast_file" or ( warn "\n!!! Can't read $blast_file...!!! " and die );

my $total = 0;

my %n50 = ();
my %n85 = ();
my %n99 = ();
my %alignements = ();
while (<BLAST>) {
    chomp;
    my ($qid, $sid, $id, $aln_len, $mis, $gaps, $q_start, $q_end, $s_start, $s_end, $eval, $bit) = split /\s*\t\s*/, $_;
    if (lc($e_val) eq "0.0") { $e_val = 0; print STDERR "eval $e_val\n"; }

    $qid =~ /([\d\.]+)\_.*/;
    my $cluster = $1;

    if (! defined($n50{$cluster})) { $n50{$cluster} = 0 }
    if (! defined($n85{$cluster})) { $n85{$cluster} = 0 }
    if (! defined($n99{$cluster})) { $n99{$cluster} = 0 }
    if ( ! defined($alignements{$cluster})) {
	$alignements{$cluster} = [ ($id, $aln_len, $gaps, $q_start, $q_end, $eval) ];
    }
    else {
	my ($c_id, $c_aln_len, $c_gaps, $c_q_start, $c_q_end, $c_eval) = @{$alignements{$cluster}} ;
	if ($eval != 0 && $eval < $c_eval) {
	    $alignements{$cluster} = [ ($id, $aln_len, $gaps, $q_start, $q_end, $eval) ];
	}
    }
    $total++;
}
print STDERR " processed $total alignements for ", scalar(keys %alignements), " clusters\n";

foreach my $cluster (keys %alignements) {
    my ($id, $aln_len, $gaps, $q_start, $q_end, $eval) = @{$alignements{$cluster}} ;
    my $q_len = abs($q_end - $q_start);
    my $t = $aln_len/$q_end; # should approximate the prot len

    if ($t > 0.5) { $n50{$cluster} = $n50{$cluster} + 1; }
    if ($t > 0.85) { $n85{$cluster} = $n85{$cluster} + 1; }
    if ($t > 0.99) { $n99{$cluster} = $n99{$cluster} + 1; }

    print STDOUT "$cluster\t", $t, "\t", $n50{$cluster}, "\t", $n85{$cluster}, "\t", $n99{$cluster}, "\n";
}
