#!/usr/local/bin/perl -w

# This program receives a tsv file containing blast with header as
# MUTATION_ID POSITION  MUTATION  RELATIVE_POSITION OVERLAPPED_GENES  AFFECTED_GENES


# gets parameters
my $fasta     = $ARGV[0]; # 
my $output    = $ARGV[1]; #


# Checks if arguments were given
if (! $ARGV[0])
{
        print "No fasta file  parameter given\n";
        exit;
}elsif (! $ARGV[1])
{
        print "No out file  parameter given\n";
        exit;
}

# Open the given files
open (IN,  "<$fasta")      or die "Major problem: cannot open $fasta for reading: $!";
open (OUT,  ">$output")     or die "Major problem: cannot open $output for writing: $!";



# starts hash
my %blastHash;
my %fastaHash;



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
  my @fields              = split /\s/, $id;
  my $transcriptId        = $fields[0];
  my $source              = $fields[1];
  my $info                = $fields[2];
  my $geneBiotype         = $fields[4];
  my $transcriptBiotype   = $fields[5];

  # extract chr and pos info from transcript
  my @infoTranscript      = split /:/, $info;
  my @auxChr              = split /\-/, $infoTranscript[2];
  my $chr                 = "chr" . $auxChr[0];
  my $start               = $infoTranscript[3];
  my $end                 = $infoTranscript[4];
  my $strand              = $infoTranscript[5];

  if($strand == 1){
    $strand = "+";
  }else{
    $strand = "-";
  }

  #chr17 HAVANA  transcript  31557126  31559110  . + . gene_id "ENSG00000230471.1"; transcript_id "ENST00000428118.1";

  # order start and end asc
  if($end < $start){
    my $startAux  = $start;
    $start        = $end;
    $end          = $startAux;
  }


  print OUT "$chr\tHAVANA\ttranscript\t$start\t$end\t.\t$strand\t.\tgene_id \"$transcriptId\"\; transcript_id \"$transcriptId\"\;\n";

}



#closes files
close IN;
close OUT;
