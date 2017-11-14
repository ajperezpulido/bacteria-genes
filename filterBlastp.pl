#!/usr/bin/perl
use strict;

my $ID = $ARGV[1] || 30;
my $CO = $ARGV[2] || 85;

open blastp, $ARGV[0];
while (<blastp>) {
  chomp;

  my ($qid, $sid, $slen, $qlen, $qcov, $pid) = split/\t/;
  print "$sid\n" if $pid >= $ID && $qcov >= $CO;
}
close blastp;

exit;
