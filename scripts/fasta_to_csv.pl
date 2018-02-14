#!/usr/local/bin/perl -w

# This program receives a tsv file containing mutations with header as
# MUTATION_ID POSITION  MUTATION  RELATIVE_POSITION OVERLAPPED_GENES  AFFECTED_GENES


# gets parameters
my $fasta   = $ARGV[0]; # 
my $output  = $ARGV[1]; # output from txCdsPredict


# Checks if arguments were given
if (! $ARGV[0])	
{
        print "No mutation file  parameter given\n";
        exit;
}elsif (! $ARGV[1])
{
        print "No output parameter given\n";
        exit;
}


# Open the given files
open (IN,   "<$fasta")  or die "Major problem: cannot open $fasta for reading: $!";
open (OUT,  ">$output") or die "Major problem: cannot open $output for writing: $!";

  print OUT "id;transcriptId;chr;start;end;sequence;geneBiotype;transcriptBiotype\n";

$counter = 0;
# reads fasta
while (<IN>) {
  # go backs one line to analyse if starts with >
  seek(IN, -length($_), 1);

  # starts sequence variable as empty
  my $sequence = "";

  # reads sequence identifier
  while(<IN>) {
    if (m/^>(.+)\n/)	
    {
      $id = $1;
      last;
    }
  }

  # reads sequence until find other identifier
  while(<IN>) {
    if (m/^>.+\n/)
    {
      last;
    }
    $sequence .= $_;
  }

  $sequence =~ s/\s//g;
  # splits ids by \t
  my @fields              = split /\s/, $id;
  my $transcriptId        = $fields[0];
  my $source              = $fields[1];
  my $info                = $fields[2];
  my @auxGene             = split /:/, $fields[4];
  my $geneBiotype         = $auxGene[1];

  my @auxBiotype             = split /:/, $fields[5];
  my $transcriptBiotype   = $auxBiotype[1];

  # extract chr and pos info from transcript
  my @infoTranscript      = split /:/, $info;
  my $chr                 = $infoTranscript[2];
  my $start               = $infoTranscript[3];
  my $end                 = $infoTranscript[4];
  my $strand              = $infoTranscript[5];

  # order start and end asc
  if($end < $start){
    my $startAux  = $start;
    $start        = $end;
    $end          = $startAux;
  }

if($strand eq '1'){
	$strand = '+';
}elsif($strand eq '-1'){
	$strand = '-';
}

  print OUT "$counter;$transcriptId;$chr;$start;$end;$strand;$sequence;$geneBiotype;$transcriptBiotype\n";

}





#closes files
close IN;
close OUT;
