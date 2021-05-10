#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Copy qw/mv cp/;
use File::Find qw/find/;


# Different logging functions to STDERR
sub logmsg   {print STDERR "@_";} # simple message
sub logmsgLn {logmsg "@_\n"; }    # newline
sub logmsgBg {logmsg "$0 @_"}     # with script name (BG: beginning of line)

exit main();

sub main{

  my $settings={};
  GetOptions($settings, qw(genome=s --sequence-dir=s help)) or die $!;
  die usage() if($$settings{help});
  $$settings{'sequence-dir'}//="library";

  mkdir "taxonomy";
  mkdir $$settings{'sequence-dir'};
  mkdir $$settings{'sequence-dir'}."/gtdb";


  my $assemblyId = $$settings{'genome'};
  $assemblyId=~s/^RS_//;   # remove prefix RS_
  $assemblyId=~s/^GB_//;   # remove prefix RS_
  $assemblyId=~s/\.\d+$//; # remove version
   # Download the genome with the last taxid
  my $filename = $$settings{'sequence-dir'}."/gtdb/$assemblyId.fna";
  logmsgBg "  finding it ($assemblyId)...";
  if(-e $filename && (stat($filename))[7] > 0){
    logmsgLn "file present, not downloading again.";
    next;
  }

  # Copy or download the file
  for my $i (1..10){
    logmsg "from NCBI using esearch ($assemblyId)...";
    my $sysout = system("esearch -db assembly -query $assemblyId  | elink -target nuccore  | efetch -format fasta  > $filename");
    my $filesize = (stat "$filename")[7];
    if($sysout == 0 && $filesize > 1000){
      last;
    }
  }
  my $filesize = (stat "$filename")[7];
  if($filesize < 1000){
    exit 1;
  }
    
 

  logmsgLn;

  return 0;
}


