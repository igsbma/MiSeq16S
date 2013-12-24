#! /usr/bin/perl

=head1 NAME

  stitching_pipeline.pl

=head1 DESCRIPTION

  Given tag clean fasta files of S1 and S2 sequences stitching of the resulting sequences.

  1. Building Multiple Sequence Alignment of all S1 and S2 sequences which have paired ends
  2. Stitching pair end reads using the above alignment


  NOTES

  1. This script assumes that the amplicon region is 319F_806R. For different
  regions change the value of $refAlignment.


  2. Consider filtering sequences based on their length before executing this script.
  Ex.
  seqLenDistr.pl -i em09_fwd_seqs.fna -o em09_fwd_seqs.seqLen
  seqLenDistr.pl -i em09_rev_seqs.fna -o em09_rev_seqs.seqLen

  seqLen_filtering.pl -m 215 -i em09_fwd_seqs.fna -o em09_fwd_seqs.fa
  seqLen_filtering.pl -m 200 -i em09_rev_seqs.fna -o em09_rev_seqs.fa


=head1 SYNOPSIS

  stitching_pipeline.pl -i <S1 file> -j <S2 file> -o <output dir> [Options]

=head1 OPTIONS


=over

=item B<--S1-file, -i>
  S1 file.

=item B<--S2-file, -j>
  S2 file.

=item B<--output-dir, -o>
  Output directory.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<--print-all, -a>
  Output fasta file contains stitching for _all_ pair ends, not only for the
  representatives from 100% identity clusters.

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  cd /usr/local/projects/pgajer/projects/ContraceptiveStudy/RAV247_252

  stitching_pipeline.pl -i S1.contra_tagclean.fasta -j S2.contra_tagclean.fasta -o contra_stitch

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Parallel::ForkManager;

$OUTPUT_AUTOFLUSH = 1;

my $MOTHUR_BIN  = "/usr/local/projects/pgajer/bin/mothur";
##my $refAlignment = "/usr/local/projects/pgajer/speciateIT/spp-data/Lactobacillus/pplacer_319F_806R_nr.filter.fasta";
my $refAlignment = "/usr/local/projects/pgajer/lib/silva_db/silva.bacteria.fasta";

my $nProc = 2;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "S1-file|i=s"    => \my $s1File,
  "S2-file|j=s"    => \my $s2File,
  "output-dir|o=s" => \my $outDir,
  "print-all|a"    => \my $printAll,
  "dry-run"        => \my $dryRun,
  "help|h!"        => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$s1File)
{
  print "ERROR: Missing S1 file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$s2File)
{
  print "ERROR: Missing S2 file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outDir)
{
  print "ERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

####################################################################
##                               MAIN
####################################################################

my $startRun = time();

my @sFiles = ($s1File, $s2File);

# creating output directory
my $cmd = "mkdir -p $outDir";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $logFile = "$outDir/stitching_pipeline.log";
open LOG, ">$logFile" or die "Cannot open $logFile for writing: $OS_ERROR\n";

print "--- Building qiime seq ID => fastq seq ID table. fastq IDs identify pair ends ... ";
my $s1IdsTbl = "$outDir/s1.idsTbl";
my $s2IdsTbl = "$outDir/s2.idsTbl";

my @sIdsTbls = ($s1IdsTbl, $s2IdsTbl);

my $pm = new Parallel::ForkManager($nProc);
foreach my $i (0..1)
{
  my $sFile = $sFiles[$i];
  my $sIdsTbl = $sIdsTbls[$i];
  $pm->start and next; # do the fork
  $cmd = "extract_fastq_header.pl -i $sFile -o $sIdsTbl";
  print "cmd=$cmd\n" if $dryRun;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  $pm->finish; # terminate the child process
}
$pm->wait_all_children; # wait for all the processes which have been forked
print "done\n";


print "--- Identifying paired end reads\n";
## print S1 and S2 seq IDs to files

print "--- Parsing $s1IdsTbl ... ";
my %s1_val2id = read2colTbl_val2id($s1IdsTbl);
print "done\n";

print "--- Parsing $s2IdsTbl ... ";
my %s2_val2id = read2colTbl_val2id($s2IdsTbl);
print "done\n";

my $n1 = keys %s1_val2id;
my $n2 = keys %s2_val2id;

print "    Number of elements in $s1IdsTbl: $n1\n";
print "    Number of elements in $s2IdsTbl: $n2\n";

my $s1GoodIds = "$outDir/s1.goodIds";
my $s2GoodIds = "$outDir/s2.goodIds";

open OUT1, ">$s1GoodIds" or die "Cannot open $s1GoodIds for writing: $OS_ERROR\n";
open OUT2, ">$s2GoodIds" or die "Cannot open $s2GoodIds for writing: $OS_ERROR\n";
my $count = 0;
my %pairEnds;
my %s1_id2val;
my %s2_id2val;

for (keys %s1_val2id)
{
  if (exists $s2_val2id{$_})
  {
    print OUT1 $s1_val2id{$_} . "\n";
    print OUT2 $s2_val2id{$_} . "\n";
    my $id = $s1_val2id{$_} . ":" . $s2_val2id{$_};
    my @a = ($s1_val2id{$_} . ":S1", $s2_val2id{$_}  . ":S2");
    $pairEnds{$id} = \@a;
    $s1_id2val{ $s1_val2id{$_} } = $_;
    $s2_id2val{ $s2_val2id{$_} } = $_;
    $count++;
  }
}
close OUT2;
close OUT1;

print "    Number of pair ends: $count\n\n";


print "--- Generating paired end reads fasta file\n";

my $s1GoodFa = "$outDir/s1.fa";
my $s2GoodFa = "$outDir/s2.fa";
my @goodFa = ($s1GoodFa, $s2GoodFa);

print "--- Selecting good seqs from $s1File ... ";
$cmd = "select_seqs.pl -s $s1GoodIds -i $s1File -o $s1GoodFa";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
print "done\n";

print "--- Selecting good seqs from $s2File ... ";
$cmd = "select_seqs.pl -s $s2GoodIds -i $s2File -o $s2GoodFa";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
print "done\n";


print "--- Appending S1/S2 labels to seq IDs ... ";
my $s1GoodFaS = "$outDir/s1S.fa";
my $s2GoodFaS = "$outDir/s2S.fa";
my @goodFaS = ($s1GoodFaS, $s2GoodFaS);

$pm = new Parallel::ForkManager($nProc);
foreach my $i (0..1)
{
  my $inFile  = $goodFa[$i];
  my $outFile = $goodFaS[$i];
  $cmd ="rm_size_put_orientation.pl -s S" . ($i+1) . " -i $inFile -o $outFile";
  print "cmd=$cmd\n" if $dryRun;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  $pm->finish; # terminate the child process
}
$pm->wait_all_children; # wait for all the processes which have been forked
print "done\n";


print "--- Concatenating  $s1GoodFa and $s2GoodFa ... ";
my $cmbFa = "$outDir/s1s2.fa";
$cmd = "rm -f $cmbFa; cat $s1GoodFaS $s2GoodFaS > $cmbFa";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
print "done\n";


print "--- Aligning  $cmbFa ... ";
my $alignFile = "$outDir/s1s2.align";
my $filterFile = "$outDir/s1s2.filter.fasta";

my @tmp;
push @tmp, "align.seqs(candidate=$cmbFa, template=$refAlignment, flip=T, processors=4)\n";
push @tmp, "filter.seqs(fasta=$alignFile)\n";

my $scriptFile = createCommandTxt(\@tmp);

$cmd = "$MOTHUR_BIN < $scriptFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
print "done\n";


print "--- Generating sequance length files ... ";
my $noGapsFile = "$outDir/s1s2.noGaps.fa";
$cmd = "rmGaps -i $filterFile -o $noGapsFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $noGapsSeqLenFile = "$outDir/s1s2.noGaps.seqLen";
$cmd = "seqLenDistr.pl -i $noGapsFile -o $noGapsSeqLenFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $cmbSeqLenFile = "$outDir/s1s2.seqLen";
$cmd = "seqLenDistr.pl -i $cmbFa -o $cmbSeqLenFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "done\n";


## Generating em09_saLenGr2.txt in an R script
print "--- Running R script to generate file of reads that lost >2 bases in alignment ... ";
my $saLenGr2File = "$outDir/saLenGr2.txt";
my $saLenGr2File2 = "$outDir/saLenGr2v2.txt";
my $Rscript = qq~

options(stringsAsFactors = FALSE)
sLenTbl <- read.table(\"$cmbSeqLenFile\")
ids <- sLenTbl[,1]
sLen <- sLenTbl[,2]
aLen <- read.table(\"$noGapsSeqLenFile\")[,2]
saLen <- sLen-aLen
write.table(cbind(ids[saLen>2]), file=\"$saLenGr2File\", row.names=F, col.names=F, quote=F)
write.table(cbind(ids[saLen>2], saLen[saLen>2]), file=\"$saLenGr2File2\", row.names=F, col.names=F, quote=F)
~;

my $rFile = "script.R";
open OUT, ">$rFile" or die "cannot write to $rFile: $!\n";
print OUT "$Rscript";
close OUT;

$cmd = "R CMD BATCH $rFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


my $outR = $rFile . "out";
open IN, "$outR" or die "Cannot open $outR for reading: $OS_ERROR\n";
my $exitStatus = 1;

foreach my $line (<IN>)
{
  if ( $line =~ /Error/ )
  {
    print "R script crashed at\n$line";
    print "check $outR for details\n";
    $exitStatus = 0;
    exit;
  }
}
close IN;

if ( $exitStatus )
{
  print "done\n";
}


open IN, "$saLenGr2File" or die "Cannot open $saLenGr2File for reading: $OS_ERROR\n";
my $nLines = 0;
for (<IN>)
{
  $nLines++;
}
close IN;
print "    Number of seq's that lost more than 2 bases: $nLines\n";


print "--- Generating table of deleted seq's ... ";
my %delSeqIDs = readArray($saLenGr2File);
print "done\n";

print "--- Removing sequences that lost more than 2 bases in the alignment\n";
my $filter2File = "$outDir/s1s2.align.fasta";
$cmd = "select_seqs.pl -e $saLenGr2File -i $filterFile -o $filter2File";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
print "done\n";


print "--- Replace dots by '-'s in the aligned seq's ... ";
my $noDotsFile = "$outDir/s1s2.align.noDots.fa";
$cmd = "rm_dot_from_align.pl -i $filter2File -o $noDotsFile";
##$cmd = "rm_dot_from_align.pl -i $filterFile -o $noDotsFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
print "done\n";

print "--- Adding S1/S2 suffixes to IDs of s1_id2val s2_id2val ... ";
my %s1_id2val_S;
for (keys %s1_id2val)
{
  $s1_id2val_S{ $_ . ":S1" } = $s1_id2val{$_};
}

my %s2_id2val_S;
for (keys %s2_id2val)
{
  $s2_id2val_S{ $_ . ":S2" } = $s2_id2val{$_};
}
print "done\n";


print "--- parsing $noDotsFile ... ";
my %idToSeq = readFasta($noDotsFile);
my %idToSeqS1;
my %idToSeqS2;

for (keys %idToSeq)
{
  if (exists $s1_id2val_S{$_})
  {
    $idToSeqS1{$_} = $idToSeq{$_};
  }

  if (exists $s2_id2val_S{$_})
  {
    $idToSeqS2{$_} = $idToSeq{$_};
  }

  if ( !exists $s1_id2val_S{$_} && !exists $s2_id2val_S{$_})
  {
    print "ERROR: $_ not found in s1_id2val and s2_id2val\n";
    print LOG "ERROR: $_ not found in s1_id2val and s2_id2val\n";
  }
}
print "done\n";

my $n = keys %idToSeq;
$n1 = keys %idToSeqS1;
$n2 = keys %idToSeqS2;

print "    Number of elements in idToSeq: $n\n";
print "    Number of elements in idToSeqS1: $n1\n";
print "    Number of elements in idToSeqS2: $n2\n";

print LOG "    Number of elements in idToSeq: $n\n";
print LOG "    Number of elements in idToSeqS1: $n1\n";
print LOG "    Number of elements in idToSeqS2: $n2\n";

print "--- Stitching pair end reads\n";
my $sFile = "$outDir/stitched.fa";
open OUT, ">$sFile" or die "Cannot open $sFile for writing: $OS_ERROR\n";

for my $id (keys %pairEnds)
{
  my ($id1, $id2) = @{$pairEnds{$id}};

  if ( exists $idToSeqS2{ $id2 } && exists $idToSeqS1{ $id1 } )
  {
    my $fwdAseq = $idToSeqS2{ $id2 };
    my $revAseq = $idToSeqS1{ $id1 };
    # print "fwdAseq:\n$fwdAseq\n";
    # print "revAseq:\n$revAseq\n";
    my @fwdChars = split "", $fwdAseq;
    my @revChars = split "", $revAseq;
    my $n1 = scalar(@fwdChars) - 1;
    my @s;
    my $j = 0;
    for my $i (0..$n1)
    {
      if ( $fwdChars[$i] eq '-' && $revChars[$i] ne '-' )
      {
	$s[$j++] = $revChars[$i];
      }
      elsif ( $fwdChars[$i] ne '-' && $revChars[$i] eq '-' )
      {
	$s[$j++] = $fwdChars[$i];
      }
      elsif ( $fwdChars[$i] ne '-' && $revChars[$i] ne '-' )
      {
	$s[$j++] = $fwdChars[$i];
      }
    }
    my $stitch = join "", @s;
    #print "stitch:\n$stitch\n"; exit;
    #print OUT ">$sid\n$stitch\n";

    print OUT ">$id\n$stitch\n";
  }
  else
  {
    if ( !exists $idToSeqS2{ $id2 } && !exists $delSeqIDs{ $id2 } )
    {
      print "ERROR: $id2 does not exist in idToSeqS2\n";
      print LOG "ERROR: $id2 does not exist in idToSeqS2\n";
    }
    elsif ( !exists $idToSeqS1{ $id1 } && !exists $delSeqIDs{ $id1 } )
    {
      print "ERROR: $id1 does not exist in idToSeqS1\n";
      print LOG "ERROR: $id1 does not exist in idToSeqS1\n";
    }
  }
}
close OUT;
print "done\n";

close LOG;

print "--- Creating Makefile ... ";
my $makefile = "Makefile";
open OUT, ">$makefile" or die "Cannot open $makefile for writing: $OS_ERROR\n";
print OUT "clean:
	mv $outDir/stitched.fa .
	rm -rf $outDir mothur.* script.R*
	echo \"Moved stitched.fa from $outDir to CWD\"
	echo \"Deleted $outDir\"
";
close OUT;
print "done\n";


print "    Number of elements in $s1IdsTbl: $n1\n";
print "    Number of elements in $s2IdsTbl: $n2\n";
print "    Number of pair ends:             $count\n";
print "    Number of seq's that lost more than 2 bases: $nLines\n\n";

print "    Stitched sequences written to $sFile\n\n";

print "    To delete output directory and move the stitched file to the current working directory\n";
print "    Run: make clean\n\n";


my $endRun = time();
my $runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  my $timeMin = int($runTime / 60);
  my $timeSec = $runTime % 60;
  print "Completed in $timeMin:$timeSec\n"
}
else
{
  print "Completed in $runTime seconds\n"
}

exit;


####################################################################
##                               SUBS
####################################################################

# read two column table and create a table that assigns elements of the second
# column to the corresponding element of the first column
sub read2colTbl_val2id
{
  my $file = shift;

  my %val2id;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $val2id{$t}  = $id;
  }
  close IN;

  return %val2id;
}

sub createCommandTxt
{
    my (@arr) = @{$_[0]};
    my $file = "$outDir/script.txt";
    open OUT, ">$file" or die "Cannot open file $file to write: $!\n";
    map { print OUT $_ } @arr;
    print OUT "quit()\n";
    return $file;
}

sub readFasta
{
  my $inFile = shift;
  my %idToSeq;
  open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
  $/ = ">";
  my $junkFirstOne = <FASTA>;
  while (<FASTA>)
  {
    chomp;
    my ($id,@seqlines) = split /\n/, $_;
    my $seq = join '', @seqlines;
    $idToSeq{$id} = $seq;
  }
  $/ = "\n";
  close FASTA;

  return %idToSeq;
}

# read table with one column
sub readArray
{
  my $file = shift;
  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    $tbl{$_}=0;
  }
  close IN;

  return %tbl;
}

exit;
