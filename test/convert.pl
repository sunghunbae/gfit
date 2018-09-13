#!/usr/bin/perl

$percent_error = 3;

if ($#ARGV != 0) {
  print "\n\tUsage: col2dat file\n\n";
  exit;
  }

open (INPUT,$ARGV[0]) or die "cannot open $ARGV[0]\n";
@file = <INPUT>;
close INPUT;
foreach $line (@file) {
  @c = split /\s+/, $line;
  for ($i= 0; $i< $#c ;$i++ ) {
    @column[$i] .= $c[$i] . " "; }
  $ncols = $#c;
  }
for($i=1; $i < $ncols; $i++) {
  @c0 = split /\s+/, @column[0];
  @c1 = split /\s+/, @column[$i];
  $nrows = $#c0;

  # header
  printf("%s %s 1\n",$c0[0],$c1[0]);

  for ($j=1; $j<$nrows; $j++) {
    if ($c0[$j] > 0) {
      printf("%.5f %.5f %.5f\n",$c0[$j],$c1[$j],$c1[$j]*$percent_error/100.);
      }
    }
  printf("&\n");
  }
