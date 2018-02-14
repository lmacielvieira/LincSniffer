#!/usr/local/bin/perl -w

# This program receives a tsv file containing blast with header as
# MUTATION_ID POSITION  MUTATION  RELATIVE_POSITION OVERLAPPED_GENES  AFFECTED_GENES


# gets parameters
my $fasta     = $ARGV[0]; #
my $ids       = $ARGV[1]; # 
my $output    = $ARGV[2]; #


# Checks if arguments were given
if (! $ARGV[0])
{
        print "No fasta file  parameter given\n";
        exit;
}elsif (! $ARGV[1])
{
        print "No id file  parameter given\n";
        exit;
}elsif (! $ARGV[2])
{
        print "No out file  parameter given\n";
        exit;
}

# Open the given files
open (IN2,  "<$fasta")      or die "Major problem: cannot open $fasta for reading: $!";
open (IN,  "<$ids")      or die "Major problem: cannot open $ids for reading: $!";
open (OUT,  ">$output")     or die "Major problem: cannot open $output for writing: $!";



# starts hash
my %idsHash;
my %fastaHash;


# read tsc with blast
while (<IN>) {
  if (m/^(.+)\n/){
      $id = $1;
      $idsHash{$id} = $1;
  }
}


# reads fasta
while (<IN2>) {
  # starts sequence variable as empty
  my $sequence = "";

  # reads sequence identifier
  if (m/^>(.+)\n/)
    {
      $id = $1;
  }


  # reads sequence until find other identifier
  while(<IN2>) {
    if (m/^>.+\n/)
    {
      # go backs one line to analyse if starts with >
      seek(IN2, -length($_), 1);
      last;
    }
    $sequence .= $_;
  }

  # splits ids by \t
  my @fields              = split /\s/, $id;
  my $transcriptId        = $fields[0];
  my $source              = $fields[1];
  my $info                = $fields[2];
  my $geneBiotype         = $fields[4];
  my $transcriptBiotype   = $fields[5];

  # extract chr and pos info from transcript
  my @infoTranscript      = split /:/, $info;
  my $chr                 = $infoTranscript[2];
  my $start               = $infoTranscript[3];
  my $end                 = $infoTranscript[4];

  # order start and end asc
  if($end < $start){
    my $startAux  = $start;
    $start        = $end;
    $end          = $startAux;
  }

  # removes space from sequence and calculates length
  $sequence =~ s/\s//g;

  if(exists $idsHash{$transcriptId}){
    print OUT ">$id\n$sequence\n";
  }
}



#closes files
close IN;
close IN2;
close OUT;

