#!/usr/bin/perl -w
# --------------------------------------------------------------------------------
#
#  Filipe Ribeiro 2009-05-26 ribeiro@broad.mit.edu
#
#  Utility to recursively parse all 'makeinfo/<target>.out.<module> files and 
#  calculate cumulative run time. 
#
#    Uses:
#     start time: from first lines of files
#     end time: from file time stamp
#
# --------------------------------------------------------------------------------

# ---------------------------
#  time conversion constants
# ---------------------------

my $d_h = 24;           my $h_d = 1.0 / $d_h;
my $h_m = 60;           my $m_h = 1.0 / $h_m;
my $m_s = 60;           my $s_m = 1.0 / $m_s;
my $h_s = $h_m * $m_s;  my $s_h = 1.0 / $h_s;
my $d_s = $d_h * $h_s;  my $s_d = 1.0 / $d_s;





print "#---- finding 'makeinfo' files\n";
my @dns = recursive_find_fn(".", "makeinfo");
map print("#     $_\n"), @dns;



my %t0 = ();
my %t1 = ();


foreach my $dn (@dns) 
{
    my @fns = `ls -at $dn`;
    if (@fns > 2) # exclude directories with '.' and '..' only 
    {
        chomp @fns;

        my $ttot = 0;
        print "#---- parsing '$dn' for times\n";
        print "#\n";
        print "#     t_start(s)   t_run(s)          t_run     module.pid\n";
        print "#\n";
        
        foreach my $fn (@fns)
        {
            if ($fn !~ /\.mm\./ and $fn =~ /\.out\./) 
            {
                my ($modn, $t0, $t1) = module_timing("$dn/$fn");
                
                if (!exists $t0{$modn}) 
                {
                    $t0{$modn} = $t0;
                    $t1{$modn} = $t1;
                    
                    printf "      %10d %10d = %12s     %-30s\n", $t0, ($t1-$t0), timef($t1-$t0), $modn;
                    #printf "      %10d %10d %10d %-30s\n", $t0, $t1, $t1-$t0, $modn;
                    
                    $ttot += $t1 - $t0;
                }
            }
        }
        print "#\n";
        printf "#     total time = %s\n", timef($ttot);
        print "#\n";
    }
}






sub timef
{
    my ($t_s) = @_;
    return "-"  if ($t_s <= 0);
    my $d = int($t_s * $s_d);
    $t_s -= $d * $d_s;
    my $h = int($t_s * $s_h);
    $t_s -= $h * $h_s;
    my $m = int($t_s * $s_m);
    
    return sprintf("%2dd %02dh %02dm", $d, $h, $m) if ($d > 0);
    return sprintf("%2dh %02dm", $h, $m) if ($h > 0);
    return sprintf("%02dm", $m);
}




sub module_timing
{
    my ($fn) = @_;
    
    my @lines = `head -3 $fn`;

    my ($date, $pid) = $lines[2] =~ /^(.+)\s+run\s+\(pid=(\d+)\)/; 
    my ($modn) = $fn =~ /out\.([^\.]+)$/;
    $modn .= ".$pid";
    
    my $t0 = `date +"\%s" --date="$date"`; chomp $t0;
    my $t1 = (stat($fn))[9];

    return ($modn, $t0, $t1);
}






sub recursive_find_fn 
{
    my ($dn, $fn) = @_;
    #print "$dn - $fn\n";
    
    my @dns = ();

    foreach my $fn0 (`ls $dn`)
    {
        chomp $fn0;
        my $ffn0 = "$dn/$fn0";

        if (-d $ffn0 and 
            $fn0 ne "seed" and 
            $fn0 ne "RawData" and 
            $fn0 !~ /^psf/)
        {
            #print "$dn ------ $fn0\n";
            push @dns, recursive_find_fn($ffn0, $fn);
        }
        
        if ($fn0 eq $fn)
        {
            push @dns, $ffn0;
        }
    }
    return @dns;
}


