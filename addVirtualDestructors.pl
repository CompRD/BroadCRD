#!/usr/bin/perl -w 

#Add virtual destructors to classes that are missing them.
#The input file should contain the warning lines from
#gcc, sorted and uniqued.


while (<>) {
  if (/\s*(\S+):\d+: warning: '(\S+)\s+(\S+)' has virtual functions but non-virtual destructor/) { 
    print "$1, $3\n";
    $fname = $1;
    $classname = $3;
    if ($classname =~ /\S+::(\S+)/) { $classname = $1; }
    addVirtualDestructor($fname, $classname);
  }
}
exit;



sub addVirtualDestructor {
  my ($fname, $classname) = @_;
  my $bcount=0;
  my $inclass=0;
  open(IN, $fname);
  open(OUT, "> $fname.fixed");

  while (<IN>) {
    if ((/\s*class\s+(\S+)/) || (/\s*struct\s+(\S+)/)) {
      if ($1 eq $classname) { 
	$inclass = 1;
	$bcount = 0;
	$dx=0;
      }
    }
    if (!$inclass) { print OUT; next; }

    if (/\~$classname/) { $dx = 1; }
    if (/{.*{/ || (/}.*}/))  { die "More than one bracket in line $_"; }
    if(/{/) { $bcount++; }
      if(/(.*)}(.*)/) {
	$bcount--;
	$before = $1; $after=$2;
	if (0 == $bcount && 0 == $dx) {
	  $inclass = 0;
	  if ($before !~ /^\s*$/) { print OUT "$before\n"; }
	  print OUT "public:\n  virtual ~$classname()\{\}\n\}$after\n"
	}
	else { print OUT; }
      }
    else { print OUT; } 
  }
  system("cp $fname.fixed $fname");
}
