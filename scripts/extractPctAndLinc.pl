#!/usr/local/bin/perl -w

# receives the fasta from ensembl
# as first parameter
my $linc = $ARGV[0];
my $pct  = $ARGV[1];

# Open the fasta output file
$saida_linc = $ARGV[2];
$saida_pct  = $ARGV[3];

# checks if the fasta file was given as parameter
if (! $ARGV[0])
{
	print "Error opening linc fasta file\n";
	exit;
}

# checks if the fasta file was given as parameter
if (! $ARGV[1])
{
	print "Error opening pct fasta file\n";
	exit;
}

# checks if the fasta file was given as parameter
if (! $ARGV[2])
{
	print "Linc Output file parameter missing\n";
	exit;
}

# checks if the fasta file was given as parameter
if (! $ARGV[3])
{
	print "Pct Output file parameter missing\n";
	exit;
}


# opens files for reading and writing
open (IN, "<$linc") or die "Major problem: cannot open $linc for reading: $!";
open (OUT, ">$saida_linc") or die "Major problem: cannot open $saida_linc for writing: $!";

$counter = 0;

# reads fasta
while (<IN>) {
    # starts sequence variable as empty
  my $sequence = "";

  # reads sequence identifier
  if (m/^>(.+)\n/)
    {
      $id = $1;
  }


  # reads sequence until find other identifier
  while(<IN>) {
    if (m/^>.+\n/)
    {
      # go backs one line to analyse if starts with >
      seek(IN, -length($_), 1);
      last;
    }
    $sequence .= $_;
  }


  # splits ids by \t
  my @fields = split /\s/, $id;
  my $source = $fields[1];
  my $gene_biotype = $fields[4];
  my $transcript_biotype = $fields[5];

  # removes space from sequence and calculates length
	$sequence =~ s/\s//g;
	$tam = length $sequence;

  if($tam >= 200 && 
  $gene_biotype eq "gene_biotype:lincRNA" &&
  $transcript_biotype eq "transcript_biotype:lincRNA"){
  				$counter += 1;
      print OUT ">" . $id . "\n";
      print OUT $sequence . "\n";
  }


}


print "TOTAL LINC $counter\n";

close IN;
close OUT;


# opens files for reading and writing
open (IN, "<$pct") or die "Major problem: cannot open $linc for reading: $!";
open (OUT, ">$saida_pct") or die "Major problem: cannot open $saida_linc for writing: $!";

# reads fasta
while (<IN>) {
     # starts sequence variable as empty
  my $sequence = "";

  # reads sequence identifier
  if (m/^>(.+)\n/)
    {
      $id = $1;
  }


  # reads sequence until find other identifier
  while(<IN>) {
    if (m/^>.+\n/)
    {
      # go backs one line to analyse if starts with >
      seek(IN, -length($_), 1);
      last;
    }
    $sequence .= $_;
  }

  # splits ids by \t
  my @fields = split /\s/, $id;
  my $source = $fields[1];
  my $gene_biotype = $fields[4];
  my $transcript_biotype = $fields[5];

  # removes space from sequence and calculates length
	$sequence =~ s/\s//g;
	$tam = length $sequence;

  if($tam >= 200 && $source eq "cdna:known" &&
  $gene_biotype eq "gene_biotype:protein_coding" &&
  $transcript_biotype eq "transcript_biotype:protein_coding"){
      print OUT ">" . $id . "\n";
      print OUT $sequence . "\n";
  }


}



close IN;
close OUT;
