#!/usr/local/bin/perl -w

# This program receives a fasata file and the output of txCdspredict fil (.cds)
# in order to extract features related to ORFs that characterizes lincRNAs or
# PCTs such as : orf length, orf proportion, orf start, orf and and txCdspredict
# score

# gets parameters
my $fasta = $ARGV[0]; # fasta file
my $cds = $ARGV[1];   # output from txCdsPredict
my $saida = $ARGV[2]; # output file containing orfs parameters
my $label = $ARGV[3]; # data label

# Checks if arguments were given
if (! $ARGV[0])
{
        print "No fasta file  parameter given\n";
        exit;
}elsif (! $ARGV[1])
{
        print "No cds file  parameter given\n";
        exit;
}elsif (! $ARGV[2])
{
        print "No output file  parameter given\n";
        exit;
}

# Open the given files
open (IN, "<$fasta") or die "Major problem: cannot open $fasta for reading: $!";
open (IN2, "<$cds") or die "Major problem: cannot open $cds for reading: $!";
open (OUT, ">$saida") or die "Major problem: cannot open $saida for writing: $!";

print OUT "id;Class;";


my @freqs =  ('GCGG', 'TTTT', 'TCG', 'AAAA', 'ACG', 'TTGT', 'TAT', 'TAC','GTT', 'GTG', 'AGT', 'CCGA', 'TACC', 'CGTG', 'CGCT', 'TACG','TTAG', 'CGTA', 'ACCG', 'CCGT', 'CGGT', 'CGAC', 'CGCA', 'GCGT','GTAG', 'CGTT', 'CGAA', 'GCGA', 'CGAT', 'TAGT');
foreach my $kmer (@freqs) {
  print OUT $kmer . ';';
}

# prints csv header
print OUT "ORF_Start;ORF_End;ORF_Score;ORF_Length;" .
"ORF_proportion;Transcript_Length\n";


# orfs hash
my %orfs_hash;

# reads cds file
  while (<IN2>) {
    if (m/^(\S+)\s+(\S+)\s+(\S+)\s+txCdsPredict\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\n/)
    {
      $id_orf         = $1;
      $orf_start      = $2;
      $orf_end        = $3;
      $orf_score      = $4;
      $orf_length     = abs ($3 - $2);
      
      $orfs_hash{$id_orf} = {"start" => $orf_start, "end" => $orf_end, 
      "score" => $orf_score, "length" => $orf_length};
    }
  }


# reads fasta
while (<IN>) {
  # starts sequence variable as empty
  my $sequence = "";

  # reads sequence identifier
    if (m/^>(\S+)(.+)*\n/)
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

  
  $transcript_length = length $sequence;

  $orf_proportion =  $orfs_hash{$id}{"length"} / $transcript_length;
  
   #PUT SEQUENCE IN UPPER CASE
   $sequence = uc "$sequence";

    print  OUT "$id;$label;";
    foreach my $kmer (@freqs) {
       my @c_gc = $sequence =~ /$kmer/g;
       my $count_gc = @c_gc;
       print OUT $count_gc . ';';
      
    }
  print  OUT "$orfs_hash{$id}{'start'};$orfs_hash{$id}{'end'};$orfs_hash{$id}{'score'};$orfs_hash{$id}{'length'};$orf_proportion;$transcript_length\n";

}


#closes files
close IN;
close IN2;
close OUT;
