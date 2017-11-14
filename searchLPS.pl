#!/usr/bin/perl
use strict;

my $O = $ARGV[1] || "i";

open db, $ARGV[0];
while (<db>) {
  chomp;
  my $h = 0;
  my @c;

  $_ = lc $_;
  #next if     $_ =~ /(transport|resistance|response to antibiotic)/; # for bactericin
  #next unless $_ =~ /(transport|resistance|response to antibiotic)/; # for resistance

  my ($id, @co) = (split/[\t]/)[0,3,5,6,7];
  foreach my $co (@co) {
    next if !$co;
    my (@terms) = split/;/, $co;
    @c = (@c, @terms);
  }
  
  foreach my $c (@c) {
    if ($c =~ /(lps|lipopolysaccharide)/) {
      my (@t) = split/[^0-9a-z]/, $c;
      foreach my $t (@t) {
        if ($t =~ /^lps$/ || $t =~ /lipopolysaccharide/) {
          #next if ($t =~ /(scin|lectin)/); # Staphylococcus virulence & Lectin
          print "$t\n" if $O eq "t"; # .terms
          $h = 1;
        }
      }
    }
  }
  $id = uc $id;
  print "$id\n" if $h == 1 && $O eq "i"; # .id
}
close db;

exit;
