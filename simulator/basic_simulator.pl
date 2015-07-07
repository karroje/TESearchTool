#!/usr/bin/perl -w
use strict;

# Generate a .mal file and a "genome" with one repeat from the model
# Model description:
#   * Creates an mal file, and a simulated genome containing one repeat at a specified position.  End bases must be present.
#   * Does NOT truncate repeats; all repats full length.
#   * Insertion characters are picked according to background rate
# Input: A parameter file.  See sim_params.txt for example.
# Output description:
#   * $ARGV[1].mal: A standard .mal file generated from the model.
#   * $ARGV[2].fa: A fast a file containing the genome.  Info. line contains the position of the repeat in the genome
#                  and the "alignment" of this repeat as it would be represented in the .mal file.
my $rng_seed;
my $cut_ends = 0;
while ($ARGV[0] =~ /-{1,2}(\w+)/) {
    my $switch = $1;
    my $full_switch = shift @ARGV;

    if (grep {$switch eq $_} qw(s seed)) {
	$rng_seed = shift @ARGV;
    }
    elsif (grep {$switch eq $_} qw(c cut)) {
	$cut_ends = 1;
    }
    else {
	die "Bad switch: $full_switch\n";
    }
}

srand($rng_seed) if $rng_seed;

my $param_file = $ARGV[0];
my $mal_file = $ARGV[1];
my $genome_file = $ARGV[2];

$mal_file = $mal_file.".mal" unless $mal_file =~ /\.mal$/;
$genome_file = $genome_file.".fa" unless $genome_file =~ /\.fa$/;

my @substitution_matrix;  # Holds the cumulative substitution rates
my @transition_matrix;    # Cumulative
my @ancestor_rates;       # Cumulative
my @background_distribution;  #Cumulative
my $ancestor_length;
my $num_modern;
my $genome_length;
my $repeat_position;

my @alph = qw(A C G T);
my @alph_lc = qw(a c g t);
my %index = (A => 0, C => 1, G => 2, T => 3);

########### Read in parameters
my $fp;
open($fp, $param_file) || die("Cannot read file: param_file");
my $temp;

# Read in substitution rates
$temp = next_line($fp);
$temp eq "substitution rates" || die("Bad param. format: substitution rates ($temp)\n");
@substitution_matrix = read_distribution_matrix($fp, 4);

# Read in tansition rates
$temp = next_line($fp);
$temp eq "transition rates" || die("Bad param. format: transition rates ($temp)\n");
@transition_matrix = read_distribution_matrix($fp, 3);

# Read in ancestor rates
$temp = next_line($fp);
$temp eq "ancestor distribution" || die("Bad param. format: ancestor rates ($temp)\n");
my @M = read_distribution_matrix($fp, 1);
@ancestor_rates = @{$M[0]};

#Read in background distibution
$temp = next_line($fp);
$temp eq "background distribution" || die("Bad param. format: background rates ($temp)\n");
@M = read_distribution_matrix($fp, 1);
@background_distribution = @{$M[0]};

#Read in other variables
$temp = next_line($fp);
$temp eq "ancestor length" || die("Bad param. format: ancestor length ($temp)\n");
$ancestor_length = next_line($fp);

$temp = next_line($fp);
$temp eq "number of modern" || die("Bad param. format: number of modern\n");
$num_modern = next_line($fp);

$temp = next_line($fp);
$temp eq "genome length" || die("Bad param. format: genome rates ($temp)\n");
$genome_length = next_line($fp);

$temp = next_line($fp);
$temp eq "repeat position" || die("Bad param. format: repeat position ($temp)\n");
$repeat_position = next_line($fp);

# Generate mal
my $mfp;
open($mfp, ">", $mal_file);
my $ancestor = generate_seq($ancestor_length, \@ancestor_rates);
print $mfp "ancestor\t$ancestor\n";
#my ($start, $stop) = (0, length($ancestor)-1);
for ((1..$num_modern)) {
    my ($descendent, $left_end, $right_end) = generate_descendant($ancestor, \@substitution_matrix, \@background_distribution, \@transition_matrix, $cut_ends);
    print $mfp join "\t", ("chr?", "?", "?", $left_end, $right_end, $descendent);
    print $mfp "\n";
}
close($mfp);

# Genorate genome
my $gfp;
open($gfp, ">", $genome_file);
my $g1 = generate_seq($repeat_position, \@background_distribution);
my ($rpt,$start,$stop) = generate_descendant($ancestor, \@substitution_matrix, \@background_distribution, \@transition_matrix, $cut_ends);
my $g3 = generate_seq($genome_length - $repeat_position, \@background_distribution);
print $gfp ">$repeat_position\t$start\t$stop $rpt\n";
$rpt =~ s/-//g;
$rpt = uc($rpt);
print $gfp "$g1$rpt$g3";




sub generate_descendant {
    my $ancestor_seq = $_[0];
    my @M = @{$_[1]};
    my @bg = @{$_[2]};
    my @transitions = @{$_[3]};
    my $cut_ends = $_[4];
    my $seq = "";

    my $current_state = 0;   # 0 = match, 1 = insert, 2 = delete
    my $current_base = 0;

    my ($left_end, $right_end) = (0, length($ancestor_seq));
    if ($cut_ends) {
	$left_end = int(rand(length($ancestor_seq)));
	$right_end = int(rand(length($ancestor_seq) - 1));

	if ($right_end >= $left_end) {
	    $right_end++;
	}
	else {
	    ($left_end, $right_end) = ($right_end, $left_end);
	}
	$ancestor_seq = substr($ancestor_seq, $left_end, $right_end - $left_end);
    }

    while (1) {
	# Pick base
	if ($current_state == 0 || $current_base == length($ancestor_seq) - 1) {
	    $seq .= $alph[choose_from_cumulative($M[$index{substr($ancestor_seq, $current_base, 1)}])];
	}
	elsif ($current_state == 2) {
	    $seq .= "-";
	}
	else {
	    $seq .= $alph_lc[choose_from_cumulative(\@bg)];
	}
	last if $current_base == length($ancestor_seq) - 1;
	   
	# Pick next state
	$current_state = choose_from_cumulative($transitions[$current_state]);
	$current_base++ if $current_state != 1;  # Move to a new base in the ancestor if you are not in an insert state
    }

    return ($seq, $left_end, $right_end-1);
}

	   


sub generate_seq {
    my $size = $_[0];
    my @dist = @{$_[1]};
    return join "", map {$alph[choose_from_cumulative(\@dist)]} (1..$size);
}

sub next_line {
    my $fp = $_[0];
    my $temp;
    do {
	chomp($temp = <$fp>);
    } while ($temp =~ /^(\s*$)|\#/);
    $temp;
}


sub read_distribution_matrix {
    my @matrix;
    my ($fp, $size) = @_;
    my $row;
    for my $i ((1..$size)) {
	my @cumulative = distribution_to_cumulative(split /\s+/, next_line($fp));
	push @matrix, \@cumulative;
    }
    return @matrix;
}


sub distribution_to_cumulative {
    my @cumulative;
    my @distribution = @_;
    my $total = 0;
    for (@distribution) {
	$total += $_;
	push(@cumulative, $total);
    }

    $cumulative[-1]==1 || die("Distribution does not sum to 1\n");
    return @cumulative;
}

sub choose_from_cumulative {
    my @D = @{$_[0]};
    my $r = rand();
    my $i = 0;
    $i++ while ($D[$i] < $r);
    return $i;
}

#DEBUGGING
sub print_matrix { 
    print join "\n", map {join "\t", @$_} @{$_[0]};
    print "\n";
}


