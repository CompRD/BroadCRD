#!/util/bin/perl -w
use strict;
use FileHandle;

# list of global variables
use vars qw($bVerbose
            $bForceBuild
	    $bForceRun
            $specifiedProject
	    $projectsSource
            $emailRecipients
            $bHelp
            $rebuildFrequency
            %finishedProjects
	    %badProjects
	    $bNoRebuild
	    $buildOptions
	    $assemblyOptions
	    $testAssemblyOptions
	    $refAssemblyOptions
            $bRepeatProjects
	    $bRequireKnownContigs
	    $bSaveAllProjects
            $workingDir
            $refBuild
            $testBuild
	    $bNewCodeStats
	    $bNewCodeStatsCumulative);

main();

sub PrintOptions
{
    print << "EOF";
Usage: regress_test.pl [options]
Options:
--help                 Display this information.
--email=<recipients>   Space delimited list of users to receive error emails.
--force_run            Override checks to prevent multiple runs in same working dir.
--force_build          Force a new check out and complete rebuild of source.
--iterations=<integer> Number of iterations to run.  Default is infinite.
--max_cpu_load=<number>
                       Specify the maximum CPU load.  The script will not
                           run if the current CPU load is greater than
                           max_cpu_load.  Default is no limit.
--project=<dir>        Specify a project to compare.
--projects_source=<dir> Specify a directory from which projects are selected.
--working_dir=<dir>    Specify a directory where builds and logs are stored.
--ref_build=<dir>      Specify the reference build for comparison.  The default
                           is <working_dir>/build/Arachne.reference.  If not
                           specified and Arachne.reference does not exist, an
                           assembly is done with the test build, but no comparison
                           is performed.
--test_build=<dir>     Specify the test build for comparison.  If none specified,
                           a new Arachne source will be checked out and built in
                           <working_dir>/build every <rebuild_frequency> minutes.
--no_rebuild           Never rebuild the test_build.
--rebuild_frequency=<minutes>
                       Number of minutes between each new build.  Default is
                           24 hours.
--build_options=<options>
                       Options to pass to 'make' command during build.
--assembly_options=<options>
                       Options to pass to the assembly command.
--test_assembly_options=<options>
                       Options to use for the test assembly.  Overrides the value 
                           of --assembly_options.
--ref_assembly_options=<options>
                       Options to use for the reference assembly.  Overrides
                           the value of --assembly_options.
--no_repeat_projects   Specify not to assemble projects that are assembled
                           during the script execution.
--save_all_projects    Ordinarily, only projects that fail the regression test are
                       copied to the projects directory.  If this option is specified,
                       copy every tested project to the projects directory.
--require_known        Require selected projects to have a known contigs file.  If
                       this is not specified, the test assembly is compared to the
                       contigs produced in the reference run.  If this is specified, 
                       both the reference and test assemblies are compared to the
                       given known contigs.
--sleep=<seconds>      Number of seconds to sleep before trying to run again
                           if the CPU load is greater than max_cpu_load.
                           Default is 5 minutes.
--verbose              Print extra run-time information.

--new_code_stats       Calculate counts for new lines of code and write results to a 
                       file in '--working_dir'
--new_code_stats_cumulative 
                       if --new_code_stats specified, produces a single file containing
		       the cumulative counts for each new line of code over all
		       iterations.

EOF
}

sub main
{
    autoflush STDOUT 1;

    # set our environment variables so we can find our executables
    #$ENV{PATH} = "/util/bin:/usr/bin";
    $ENV{ARACHNE_PRE} = "/wga/data01/WGAdata";

    use Getopt::Long;

    $bVerbose = 0;
    $bForceBuild = 0;
    $bForceRun = 0;
    $specifiedProject = "";
    $emailRecipients = "";
    $bHelp = 0;
    %finishedProjects = ();
    %badProjects = ();
    $bRepeatProjects = 1;
    $bRequireKnownContigs = 0;
    $bNoRebuild = 0;
    $bSaveAllProjects = 0;
    $rebuildFrequency = 1440;
    $workingDir = "/wga/devel/regress_test";
    $projectsSource = "/wga/scratch10/HumanProjects/Projects";
    $assemblyOptions = "";
    $testAssemblyOptions = "";
    $refAssemblyOptions = "";
    $bNewCodeStats = 0;
    $bNewCodeStatsCumulative = 0;

    my $maxCPULoad = 0;
    my $sleepDuration = 300;
    my $iterations = 0;
    my $bNoRepeat = 0;

    &GetOptions("force_build"         => \$bForceBuild,
		"force_run"           => \$bForceRun,
                "project=s"           => \$specifiedProject,
		"projects_source=s"   => \$projectsSource,
                "verbose"             => \$bVerbose,
                "email=s"             => \$emailRecipients,
                "help"                => \$bHelp,
                "max_cpu_load=s"      => \$maxCPULoad,
                "sleep=s"             => \$sleepDuration,
                "iterations=i"        => \$iterations,
                "no_repeat_projects"  => \$bNoRepeat,
		"require_known"       => \$bRequireKnownContigs,
                "working_dir=s"       => \$workingDir,
		"no_rebuild"          => \$bNoRebuild,
		"save_all_projects"   => \$bSaveAllProjects,
                "rebuild_frequency=i" => \$rebuildFrequency,
		"build_options=s"     => \$buildOptions,
		"assembly_options=s"  => \$assemblyOptions,
		"test_assembly_options=s"  => \$testAssemblyOptions,
		"ref_assembly_options=s"   => \$refAssemblyOptions,
                "ref_build=s"         => \$refBuild,
                "test_build=s"        => \$testBuild,
	        "new_code_stats"      => \$bNewCodeStats,
	        "new_code_stats_cumulative"      => \$bNewCodeStatsCumulative);

    $ENV{PATH} .= ":$workingDir/bin";

    if ($bHelp)
    {
        PrintOptions();
        exit;
    }

    die "The working directory $workingDir does not exist.\n" 
      if ( ! -d "$workingDir" );

    foreach my $subdir ( ( 'build', 'tmp', 'projects', 'bin' ) ) 
    {
	mkdir( "$workingDir/$subdir", 0777 ) 
	  if ( ! -d "$workingDir/$subdir" );
    }

    foreach my $link ( ( 'dtds', 'e_coli', 'e_coli_transposons', 'vector' ) )
    {
	symlink( $ENV{'ARACHNE_PRE'} . "/$link", "$workingDir/$link" )
	  if ( ! -e "$workingDir/$link" );
    }

    symlink("/home/strontium/jbutler/bin/addr2line",
            "$workingDir/bin/addr2line");

    if ( $bNoRebuild && $bForceBuild )
    {
	die "You cannot specify both --force_build and --no_rebuild.\n";
    }

    if ($bNoRepeat)
    {
        $bRepeatProjects = 0;
    }

    if ($specifiedProject ne '')
    {
	$iterations = 1;
    }

    if ($bNewCodeStats && $testBuild eq "")
    {
      die "You must specify --test_build if you want to use --new_code_stats\n";
    }

    if ($bNewCodeStats && (not $bNoRebuild))
    {
      die "You must specify --no_rebuild if you want to use --new_code_stats\n";
    }

    if ($bNewCodeStats && $bNewCodeStatsCumulative)
    {
      die "At this point you can't specify both --new_code_stats and --new_code_stats_cumulative\n";
    }
    
    $testAssemblyOptions = $assemblyOptions if ( $testAssemblyOptions eq '' );
    $refAssemblyOptions = $assemblyOptions if ( $refAssemblyOptions eq '' );

    if ( ( not $bRequireKnownContigs ) &&
         $testAssemblyOptions !~ /known_contigs=/ && 
	 $refBuild ne '' )
    {
	$testAssemblyOptions .= " known_contigs=ref_contigs.fasta";
    }

    # hack to disallow multiple instances to run at the same time.
    # A more elegant solution would wrap the compilation part of the
    # script with a file lock, but there's no simple way of doing it.
    # So, this tests to see if any of the files in the working dir has
    # been modified in the last x number of seconds.
    if ( ! $bForceRun &&
	 ( time() - LastLogFileModifiedTime($workingDir) < 90 ) )
    {
        PrintMessage("Another instance is running.  Stopping.\n");
        exit;
    }

    $rebuildFrequency = $rebuildFrequency * 60;

    #  if statistics on new code desired, set up
    if ($bNewCodeStats || $bNewCodeStatsCumulative) 
    {
      #  cvs diff test build with repository - output to file in $workingDir/tmp
      PrintMessage( "Getting diffs from CVS." );
      CVSDiffsToFile($workingDir, $testBuild);
      
      #  array of modified files in $testBuildDir
      my @modified_files = ModifiedFilesList($workingDir, $testBuild);
      
      if ( ($#modified_files) < 0 ) {
	die "\nNo files in $testBuild have been modified.\n";
      }
      
      #  dependants of the modified files - output to file in $workingDir/tmp
      PrintMessage( "Finding dependants of modified files.\n" );

      my $has_dependants = FindDependants($workingDir, $testBuild, @modified_files);
      if ( not $has_dependants )
      {
	die "\nWas unable to find dependants of modified files.\n";
      }

      #  build the dependants of the modified files in $testBuild with GCOV=yes OPTIM=none
      PrintMessage( "Rebuilding dependants with gcov enabled.\n" );
      
      my $is_built = BuildWithGCOV($workingDir, $testBuild);
      if ( not $is_built  )
      {
	die "\nWas unable to rebuild with gcov.\n";
      }
    }


    my $numRuns = 0;
    while (1)
    {
        if (StopScript())
        {
            PrintMessage("Found stop file.  Quitting.\n");
            last;
        }

        my $cpuLoad = GetCPULoad();

        if ( ( $cpuLoad <= $maxCPULoad ) || ( $maxCPULoad == 0.0 ) )
        {
            PrintMessage("Starting iteration $numRuns\n");

            RunScript($workingDir,$specifiedProject);

	    ++$numRuns;
        }

        if ( ( $numRuns >= $iterations ) && ( $iterations > 0 ) )
        {
            PrintMessage("Completed $iterations iterations\n");
            last;
        }

        if ( ( $cpuLoad > $maxCPULoad ) && ( $maxCPULoad > 0.0 ) )
        {
            PrintMessage("Current CPU load is $cpuLoad, which is higher ");
            PrintMessage("than $maxCPULoad.  ");

	    Sleep( $sleepDuration );

            # check if a stop file is created while we were sleeping
            if (StopScript())
            {
                PrintMessage("Found stop file.  Quitting.\n");
                last;
            }
        }
    }

    if ($bNewCodeStatsCumulative)
    { 
      # get list of files on which to run gcov
      my @gcoved_files = ReadGCOVBuiltFile($workingDir);
      
      #  go ahead, run gcov
      PrintMessage( "Running gcov." );
      my $good_run = RunGCOV($testBuild, @gcoved_files);
      if ( not $good_run )
      {
	die "\nWas unable to run gcov.\n";
      }
      

      #  output one file of counts for each modified file of code and 
      #  each project
      PrintMessage( "Producing counts files.\n" );
      my $project = "tmp/cumulative";
      ProduceCountsFiles($workingDir, $testBuild, $project);
    }

    if ($bNewCodeStats || $bNewCodeStatsCumulative)
    {
      # remove object files for execs built with gcov
      PrintMessage( "Cleaning up.\n" );
      CleanGCOV($workingDir, $testBuild);
    }

    PrintMessage("All done!\n");
  }

sub RunScript
{
    my ($workingDir,$selectedProject) = @_;

    # get the PID of our script.  We use this to create unique directories
    # and runs.
    my $pid = $$;
    my $bStatus = 0;

    opendir(DIR, "$workingDir/build");
    my @arachneBuilds = grep { /^Arachne\./ } readdir(DIR);
    closedir(DIR);

    # get the current time
    my $timeNow = time();
    my $latestTime = 0;

    # go through all the builds and find the latest one
    my $build;
    my $testBuildDir = "";
    my $referenceBuildDir = "";

    if ( $refBuild ne "" )
    {
        PrintMessage( "Using user-specified reference build dir '$refBuild'.\n" );
        $referenceBuildDir = $refBuild;
    }
    else
    {
        foreach $build (@arachneBuilds)
	{
	    # there are two build directories: "Arachne.reference" and
	    # "Arachne.<date>".  The <date> component is in the form
	    # YYYYMMDDHHmm.  This convention serves as a way to give
	    # unique names to the build directories, as well as a
	    # timestamp.  We need to extract the date part to determine
	    # the creation time of the directory.
	    
	    # the reference directory is never changed.  the <date>
	    # directory is periodically updated and rebuilt.
	    
	    PrintMessage( "Found build $build.\n" );
	  
	    if ( $build =~ /\.reference$/ )
	    {
	        if ($refBuild eq "")
		{
		    $referenceBuildDir = "$workingDir/build/$build";
		}
	    }
	    else
	    {
	        my @fields = split /\./, $build;
		my $myTime =  StringToTime($fields[1]);
	      
		if ($myTime <= $latestTime)
		{
		    DeleteDirectory("$workingDir/build/$build");
		}
		else
		{
		    if ( $latestTime > 0 ) 
		    {
		      my $tmp = TimeToString($latestTime);
		      DeleteDirectory("$workingDir/build/Arachne.$tmp");
		    }
		    $latestTime = $myTime;
		}
	    }
	}
    }

    my $projectPre = $workingDir;
    my $projectData = "tmp";
    my $project = "";
    my $newRun = "new_run";
    my $refRun = "ref_run";
    my $logFile = "";
    my $errorMsg = "";
    my $relocDir = "";

    if ($selectedProject eq "")
    {
        ($selectedProject) = SelectProject("$projectsSource");
    }

    if ($selectedProject eq "")
    {
        PrintMessage("No more unassembled projects left.\n");
        exit;
    }

    PrintMessage( "Selected project $selectedProject.\n" );

    $project = "$projectData/$selectedProject.$$";
    $relocDir = "$projectPre/$project";

    # cheezy implementation of the C "do {} while(false)" loop
    while (1)
    {
        if ( $testBuild ne "" )
        {
            PrintMessage( "Using user-specified test build dir '$testBuild'.\n" );
            $testBuildDir = $testBuild;
        }
        # build made at least 24 hours ago or the user wants to do a build
        # using the latest source
        elsif ( ( not $bNoRebuild ) && 
		( $timeNow - $latestTime >= $rebuildFrequency || $bForceBuild ) )
        {
            # turn off force build.  Otherwise, we'll do a build on every
            # iteration!
            $bForceBuild = 0;

            if ( ( $latestTime > 0 ) &&
		 ( $timeNow - $latestTime >= $rebuildFrequency ) )
            {
                my $minutes = $rebuildFrequency / 60;
                PrintMessage("Latest build is older than $minutes minutes.  ");
            }
            PrintMessage("Creating new build...\n");

            # delete the older build
            if ( $latestTime > 0 )
            {
                my $tmp = TimeToString($latestTime);
                DeleteDirectory("$workingDir/build/Arachne.$tmp");
            }

            $logFile = "$workingDir/tmp/cvs.$$";
            my $str = TimeToString($timeNow);
            $testBuildDir = "$workingDir/build/Arachne.$str";

            ExecuteCommand("mkdir $testBuildDir", ".", "");
            if (CheckOutSource($workingDir, $testBuildDir, $logFile))
            {
                DeleteFile($logFile);
            }
            else
            {
                $errorMsg = "Failed to check out source code\n";
                last;
            }

            $logFile = "$workingDir/tmp/make.$$";
            if (BuildAssembler("$testBuildDir", $logFile))
            {
                DeleteFile($logFile);
            }
            else
            {
                $errorMsg = "Failed to build Arachne\n";
                last;
            }
        }
        else
        {
            my $buildTime = TimeToFormatedString($latestTime);

            PrintMessage("Using build made at $buildTime\n");

            my $str = TimeToString($latestTime);
            $testBuildDir = "$workingDir/build/Arachne.$str";
        }

        PrintMessage("selected project = $selectedProject\n");

        my $dir = GetDirName(StripExtension("$relocDir"));
        CopyDirectory("$projectsSource/$selectedProject",
                      $relocDir);

        if ($referenceBuildDir ne "")
        {
            # run with the older build so we can compare the two assemblies
            $logFile = "$relocDir/ref_run.log";
            if (RunAssembler($referenceBuildDir, $projectPre, $project,
                             $refAssemblyOptions, $refRun, $logFile)
                && ConvertMergedContigs($testBuildDir, $projectPre,
                                        $project, $refRun, $logFile))
            {
                DeleteFile($logFile);
            }
            else
            {
                $errorMsg = "Failed to assemble $selectedProject using ";
                $errorMsg .= " reference build\n";
                last;
            }
        }
        else
        {
            PrintMessage("There is only 1 build.  No comparison is done.\n");
        }

        # run with the newer build
        $logFile = "$relocDir/new_run.log";
        if (RunAssembler($testBuildDir, $projectPre, $project,
                         $testAssemblyOptions,
                         $newRun, $logFile))
        {
            DeleteFile($logFile);
        }
        else
        {
            my $buildTime = TimeToFormatedString($latestTime);
            $errorMsg = "Failed to assemble $selectedProject using ";
            $errorMsg .= " build checked out at $buildTime\n";
            last;
        }

        if ($referenceBuildDir ne "")
        {
            # compare the two runs
            $logFile = "$relocDir/compare.log";
            if (CompareResults($testBuildDir, $projectPre, $project,
                               $newRun, $refRun, $logFile))
            {
                DeleteFile($logFile);
            }
            else
            {
                $errorMsg = "DiffAssemblies completed with errors ";
                $errorMsg .= "on $selectedProject.";
                last;
            }
        }

       $bStatus = 1;  # set the result flag to true to indicate everything
                       # completed with no errors

        last;  # must be the last line in the while() loop!
      }

    if (!$bStatus)  # something failed
    {
        if ($emailRecipients ne "")
        {
            EmailErrors($emailRecipients, $timeNow, $selectedProject,
                        "$projectPre/$project", $errorMsg, $logFile);
        }
        else
        {
            my $strTimeNow = TimeToFormatedString($timeNow);
            PrintMessage("regression test started at $strTimeNow\n");
            PrintMessage("--------------------------------------------------");
            PrintMessage("-----------------------------\n");
            PrintMessage("Error Report:\n");
            system("cat $logFile");
        }
    }

    if ( !$bStatus || $bSaveAllProjects )
    {
        # rename the log file with a ".log" extension instead of the pid.
        my @tmp = split /\./, $logFile;
        @tmp = split /\//, $tmp[0];
        my $count = @tmp;
	if ( -e $logFile )
	{
            CopyFile($logFile, "$projectPre/$project/$tmp[$count - 1].log");
	    DeleteFile($logFile);
	}
        MoveDirectory("$projectPre/$project", "$projectPre/projects");
    }
    else
    {
        DeleteDirectory("$projectPre/$project");
	
        PrintMessage("regression test on $selectedProject completed without errors\n");
    }

    if ($bNewCodeStats)
    { 
      # get list of files on which to run gcov
      my @gcoved_files = ReadGCOVBuiltFile($workingDir);
      
      #  go ahead, run gcov
      PrintMessage( "Running gcov." );
      my $good_run = RunGCOV($testBuild, @gcoved_files);
      if ( not $good_run )
      {
	die "\nWas unable to run gcov.\n";
      }


      #  output one file of counts for each modified file of code and 
      #  each project
      PrintMessage( "Producing counts files.\n" );
      ProduceCountsFiles($workingDir, $testBuild, $project);

      # remove .da files so that each iteration has an accurate count
      my $rm_cmmd = "rm $testBuild/*.da";
      system("$rm_cmmd");     
    }
  }
  

# generic function to run a command via system()
sub ExecuteCommand
{
    my ($command, $workingDir, $logFile) = @_;
    my $oldWorkingDir = `pwd`;

    PrintMessage("$command...");
    chdir($workingDir) or die "Can't cd into $workingDir";

    # output to a log file if a file name is specified
    if ($logFile ne "")
    {
        $command = "$command > $logFile 2>&1";
    }
    my $result = system("$command");

    chdir($oldWorkingDir);
    PrintMessage("  Done!\n");

    if ($result != 0)
    {
      system("echo \"\\\"$command\\\" terminated with return code = $result\" >> $logFile");
    }

    return $result == 0;
}

sub CheckOutSource
{
    my ($workingDir, $buildDir, $logFile) = @_;
    my $str = "/util/bin/cvs -d /wga/devel/ArachneCVS co -d $buildDir Arachne";
    return ExecuteCommand($str, $workingDir, $logFile);
}

sub BuildAssembler
{
    my ($workingDir, $logFile) = @_;
    return ExecuteCommand("/util/bin/make $buildOptions assembler Fasta DiffAssemblies",
                          $workingDir, $logFile);
}

sub CreateWorkingDir
{
    my ($workingDir) = @_;

    my $bDirExists = 0;
    -d $workingDir and $bDirExists = 1;

    if (!$bDirExists)
    {
        PrintMessage("$workingDir does not exist.  I'm creating it.\n");
        system("mkdir $workingDir");
        system("chmod ug+rwx $workingDir");
        system("chmod o+rx-w $workingDir");
    }

    return 1;
}

sub RunAssembler
{
    my ($workingDir, $pre, $data, $options, $run, $logFile) = @_;

    PrintMessage( "Running assembly from '$workingDir'.\n" );
    my $command = "./Assemble PRE=$pre DATA=$data RUN=$run $options";
    return ExecuteCommand($command, $workingDir, $logFile);
}

sub ConvertMergedContigs
{
    my ($workingDir, $pre, $data, $run, $logFile) = @_;

    my $bResult = ExecuteCommand("./Fasta PRE=$pre DATA=$data RUN=$run "
                                 . "HEAD=work/mergedcontigs",
                                 $workingDir, $logFile);

    if ($bResult == 0)
    {
        return 0;
    }

    return CopyFile("$pre/$data/$run/work/mergedcontigs.fasta",
                    "$pre/$data/ref_contigs.fasta");
}

sub CompareResults
{
    my ($workingDir, $pre, $project, $newRun, $refRun, $logFile) = @_;

    my $runDir;

    foreach $runDir ( $newRun, $refRun ) {
      my $bResult = ExecuteCommand( "./RetrofitAssembly PRE=$pre DATA=$project "
				    . "RUN=$runDir/work", 
				    $workingDir, $logFile);

      return $bResult if ( $bResult != 0 );
    }

    my $command = "./DiffAssemblies PRE=$pre DATA=$project NEW_RUN=$newRun OLD_RUN=$refRun";

    if ( $bRequireKnownContigs )
    {
        $command .= " CHECK_QUALITY=False";
    }
    else
    {
        $command .= " CHECK_QUALITY=True";
    }
	
    return ExecuteCommand( $command, $workingDir, $logFile);

}

# split the date string into its components and then convert it into a time_t
sub StringToTime
{
    use Time::Local;
    my $input = $_[0];

    my $pos = 0;
    my $year = "";
    my $month = "";
    my $day = "";
    my $hour = "";
    my $minute = "";

    my $byte;
    foreach $byte (split //, $input)
    {
        if ($pos >= 10)  # minute
        {
            $minute = "$minute$byte";
        }
        elsif ($pos >= 8)  # hour
        {
            $hour = "$hour$byte";
        }
        elsif ($pos >= 6)  # day
        {
            $day = "$day$byte";
        }
        elsif ($pos >= 4) # month
        {
            $month = "$month$byte";
        }
        else  # year
        {
            $year = "$year$byte";
        }

        ++$pos;
    }

    my $time = timelocal(0, $minute, $hour, $day, $month - 1, $year - 1900);
    return $time;
}

sub TimeToString
{
    use Time::Local;
    my $input = $_[0];

    my ($sec, $min, $hour, $day, $month, $year, $wday, $yday, $isdst)
        = localtime($input);

    $year += 1900;
    ++$month;

    # put 0's in front of values with only 1 digit
    if ($month < 10)
    {
        $month = "0$month";
    }
    if ($day < 10)
    {
        $day = "0$day";
    }
    if ($hour < 10)
    {
        $hour = "0$hour";
    }
    if ($min < 10)
    {
        $min = "0$min";
    }
    my $strTime = "$year$month$day$hour$min";
    return $strTime;
}

sub TimeToFormatedString
{
    my $input = $_[0];

    my @array = split(//, TimeToString($_[0]));

    my $i = 0;
    my $len = @array;
    my $strTime = "";

    while ($i < $len)
    {
        $strTime = "$strTime$array[$i]";

        if ($i == 3 || $i == 5)
        {
            $strTime = "$strTime-";
        }
        elsif ($i == 7)
        {
            $strTime = "$strTime @ ";
        }
        elsif ($i == 9)
        {
            $strTime = "$strTime:";
        }

        ++$i;
    }
    return $strTime;    
}

sub DeleteDirectory
{
    use File::Path;

    if ($_[0] ne "")
    {
        PrintMessage("Deleting directory $_[0]\n");
        rmtree($_[0]);
    }
}

sub DeleteFile
{
    my ($file) = @_;

    if ($file ne "")
    {
        if ($bVerbose)
        {
            PrintMessage("Deleting file $file\n");
        }
        unlink($file) or warn "Can't delete file $file";
    }
}

sub SelectProject
{
    my ($projectsDir) = @_;

    opendir(DIR, $projectsDir);
    # get all entries in the directory except for . and ..
    my ( @projects ) = grep !/^\.\.?$/, readdir DIR;
    closedir(DIR);

    my $numProjects = scalar @projects;
    my $numFinishedProjects = scalar keys %finishedProjects;

    PrintMessage( "Found $numProjects projects.\n" );

    # no more unpicked projects
    if ($numFinishedProjects == $numProjects)
    {
        return "";
    }

    # select a project
    while (1)
    {
        my $selectedProject = $projects[ int rand $numProjects ];

	if ( IsBadProject( $projectsDir, $selectedProject ) ) {
	  PrintMessage( "Selected project $selectedProject, but it was bad.\n" );
	  next;
	}

        if ($bRepeatProjects == 1)
        {
            return $selectedProject;
        }

        # check if we have picked it before
        unless ($finishedProjects{$selectedProject})
        {
            # we haven't seen $selectedProject before
            $finishedProjects{$selectedProject} = 1;
            return $selectedProject;
        }
    }
    return "";
}

sub IsBadProject
{
    my ($projectsDir, $project) = @_;

    if ( $badProjects{$project} )
    {
        $badProjects{$project} = 1;
        return 1;
    }

    open( PHASE, "< $projectsDir/$project/phase" );
    my $phase = <PHASE>;
    chomp $phase;
	
    if ( $phase == 0 || $phase == -1 )
    {
	$badProjects{$project} = 1;
	return 1;
    }

    open( STATUS, "< $projectsDir/$project/run/work/stamps/status" );
    my $status = <STATUS>;
    close STATUS;
    chomp $status;
	
    if ( $status ne 'done' )
    {
	$badProjects{$project} = 1;
	return 1;
    }

    if ( $bRequireKnownContigs and
	 not -f "$projectsDir/$project/contigs.fasta" )
    {
	$badProjects{$project} = 1;
	return 1;
    }

    return 0;
}
    

sub EmailErrors
{
    my ($emailRecipients, $timeNow, $selectedProject, $projectDir,
        $errorMsg, $logFile) = @_;

    if ($emailRecipients eq "")
    {
        return 0;
    }

    my $user;
    foreach $user (split(" ", $emailRecipients))
    {
        open(MAIL, "|/usr/sbin/sendmail -oi -t -odq")
            or die "Can't open sendmail: $!\n";

        PrintMessage("Sending email to $user\n");

        print MAIL "To: $user\n";
        print MAIL "From: eyee\n";
        print MAIL "Subject: Regression test failed for project "
            . "$selectedProject\n";

        my $strTimeNow = TimeToFormatedString($timeNow);
        print MAIL "regression test started at $strTimeNow\n";
        print MAIL "Project is located at $projectDir\n\n";
        print MAIL "Error Message: $errorMsg\n\n";

        if ($logFile ne "")
        {
            print MAIL "*****************************************************";
            print MAIL "**************************\n";
            print MAIL "Error Report:\n";

            open(logFile, "< $logFile") or die "Can't open file $logFile";
            my $line = "";
            while ($line = <logFile>)
            {
                print MAIL $line;
            }
            close(logFile);
        }

        close(MAIL);
    }
}

sub PrintMessage
{
    if ($bVerbose)
    {
        print("$_[0]");
    }
}

sub GetCPULoad
{
    my $load = `uptime`;
    my @tokens = split(" ", $load);

    # remove the trailing ','
    @tokens = split(",", $tokens[9]);
    return $tokens[0];
}

# test for the existence of a file with the name "stop_regression_test".  Returns 1 if it
# exists
sub StopScript
{
    my $strStopFileName = "$workingDir/stop_regression_test";

    my $bExists = ( -f $strStopFileName );
    return $bExists;
}

sub LastLogFileModifiedTime
{
    my ($workingDir) = @_;

    $workingDir .= '/tmp';
    PrintMessage("Working dir = $workingDir\n");

    opendir(DIR, $workingDir);
    # get all entries in the directory except for . and ..
    my @dirEntries = grep !/^\.\.?$/, readdir DIR;
    closedir(DIR);

    my $latestTime = 0;
    my $entry;
    foreach $entry (@dirEntries)
    {
        my ($writeTime) = (stat("$workingDir/$entry"))[9];

        if ($latestTime < $writeTime)
        {
            $latestTime = $writeTime;
        }
    }

    my $str = ( $latestTime == 0 ? "the epoch" : TimeToFormatedString($latestTime) );
    PrintMessage("Latest file modification time is $str.\n");

    return $latestTime;
}

sub MoveDirectory
{
    my ($srcDir, $dstDir) = @_;

    ExecuteCommand("mkdir -p $dstDir", ".", "");
    ExecuteCommand("mv -f $srcDir $dstDir", ".", "");
}

sub CopyDirectory
{
    my ($srcDir, $dstDir) = @_;

    ExecuteCommand("cp -r -p $srcDir $dstDir", ".", "");
}

sub CopyFile
{
    my ($strFile, $destination) = @_;

    return ExecuteCommand("cp -f $strFile $destination", ".", "");
}

sub StripExtension
{
    my ($str) = @_;

    my @tmp = split /\./, $str;
    return $tmp[0];
}

sub GetDirName
{
    my ($dirPath) = @_;

    my @tmp = split /\//, $dirPath;
    my $count = @tmp;

    if ($count > 0)
    {
        return $tmp[$count - 1];
    }
    return "";
}

sub RenameFile
{
    my ($oldName, $newName) = @_;
    return ExecuteCommand("mv $oldName $newName", ".", "");
}

sub Sleep
{
    my ( $sleepDuration ) = @_;

    open( SLEEPLOG, ">$workingDir/tmp/sleep.$$" );

    my $wakeupTime = TimeToFormatedString(time() + $sleepDuration);
    PrintMessage("Sleeping until $wakeupTime before trying again...\n");
    print SLEEPLOG "Sleeping until $wakeupTime before trying again...\n";
    my $timeSlept = 0;
    while ( $timeSlept < $sleepDuration )
    {
	my $timeToSleep = $sleepDuration - $timeSlept;
	$timeToSleep = 60 if ( $timeToSleep > 60 );
	sleep( $timeToSleep );
	$timeSlept += $timeToSleep;
	print SLEEPLOG "Slept $timeSlept seconds...\n";
	
	if (StopScript())
	{
	    $timeSlept = $sleepDuration;
	}
    }

    close SLEEPLOG;

    unlink( "$workingDir/tmp/sleep.$$" );
}

# cvs diff a build directory with the repository, output to a file in working_dir/tmp.
sub CVSDiffsToFile
{
  my ($working_dir, $build_dir) = @_;

  chdir ( $build_dir ) or die "Can't cd into $build_dir";

  my $cmmd_str = "/util/bin/cvs -d /wga/devel/ArachneCVS -q diff > $working_dir/tmp/cvs.diffs";

  my $res = system( $cmmd_str );
  if ( $res != 0 )
  {
    PrintMessage("\n$cmmd_str terminated with return code $res \n");
  }
  return $res == 0;

}

#  open the cvs.diffs file for test_build and get the files
sub ModifiedFilesList
{
  my ($working_dir) = $_[0];

  open( cvs_diffs, "<$working_dir/tmp/cvs.diffs" ) 
    or die "Can't open $working_dir/tmp/cvs.diffs";
  
  my @found_files;
  while ( <cvs_diffs> )
  {
    chomp;
    if ( $_ =~ m/^Index:/)
    {
      s/^Index: //;
      push @found_files, $_;
    }
  }
  close(cvs_diffs);

  return @found_files;
}

    
#  find which files depend on the files that have been modified 
sub FindDependants
{
  my ($working_dir, $build_dir, @modified_files) = @_;
 
  my $cmmd_str = "make FindDependencies";
  
  my $res = system( $cmmd_str );
  if ( $res != 0 )
  {  
    PrintMessage("\n$cmmd_str terminated with return code $res \n");
    return $res == 0;
  }
 
  $cmmd_str = join " ", "$build_dir/FindDependencies -g ", @modified_files;
  $cmmd_str = join " ", $cmmd_str, "> $working_dir/tmp/dependants.all";

  $res = system( $cmmd_str );
  if ( $res != 0 )
  {  
    PrintMessage("\n$cmmd_str terminated with return code $res \n");
  }
  return $res == 0;

}


#  build the execs required with GCOV=yes OPTIM=none
#  after this ready to run the regress test
sub BuildWithGCOV
{
  my ($working_dir, $build_dir) = @_;

  chdir ($build_dir) or die "Can't cd into $build_dir";

  # clean out gcov files
  print "\n";
  system("make gcov_clean");

  my @build_list = SelectFromDependants($working_dir, $build_dir);

  #  insure recompilation
  my @clean_list;
  my @exe_list;
  my $fil;
  foreach $fil (@build_list)
  {
    #  no obj files in subdirs
    my $obj_fil;
    my $filebase = RemovePath($fil);
   
    $obj_fil = join "", $fil, ".o";
    
    push @clean_list, $obj_fil; 
    push @exe_list, $filebase;
  }

  open(gcov_built, "> $working_dir/tmp/built_with_gcov.files")
    or die "Can't open $working_dir/tmp/built_with_gcov.files";
  foreach $fil (@build_list)
  {
    print gcov_built "$fil\n";
  }
  close(gcov_built);

  my $res = system( "rm @clean_list" );
  if ( $res != 0 )
  {  
    PrintMessage("\nrm failed with return code $res \n");
  }

  #  recompile with gcov and no optimization (each executable separately in case
  #  of failure)
  my $load = GetCPULoad();  
  my $ncpu_to_use = 1;
  if ( $load < 3 )
  {
    $ncpu_to_use = 2;
  }
  my $cmmd_str = "make -j$ncpu_to_use OPTIM=none GCOV=yes";
  $cmmd_str = join " ", $cmmd_str, @exe_list;
  
  print "\n";
  print "$cmmd_str\n";

  $res = system( $cmmd_str );
  if ( $res != 0 )
  {  
     PrintMessage("\ngcov build terminated with return code $res \n");
  }

  return $res == 0;

}
  
#  run regress test before running this:
#  get the files which were built with gcov so we can run gcov on them to get the 
#  line counts
sub ReadGCOVBuiltFile
{

  my $working_dir = $_[0];
  my @built_list;
  open( gcov_files, "<$working_dir/tmp/built_with_gcov.files" ) 
    or die "Can't open $working_dir/tmp/built_with_gcov.files";
  while( <gcov_files> )
  {
    chomp;
    push @built_list, $_;
  }
  close(gcov_files);

  return @built_list;
}

# 
#  for all files compiled with GCOV=yes in the test_dir, run gcov 
#  to get line counts
sub RunGCOV
{
  my ($test_build, @gcoved_list) = @_;
  
  chdir( $test_build ) or die "Can't cd into $test_build";

  my $fil;
  my $res;
  foreach $fil (@gcoved_list)
  {
    my $no_path = RemovePath($fil);
    my $cmmd_str = "gcov $no_path";

    print "\n";
    print "$cmmd_str\n";

    $res = system( $cmmd_str );
    if ( $res != 0 )
    { 
      PrintMessage("\n$cmmd_str terminated with return code $res \n");
    }
  }

  return $res == 0;
}

sub ProduceCountsFiles
{
  my ($working_dir, $test_build, $project) = @_;

  #  filenames from cvs.diffs and modified-line ranges
  my %cvs_info = ParseCVSDiffFile($working_dir, $test_build);

  #  strip "tmp/" off the project name
  my ($prfx, $tastyproject) = split(/\//, $project);  

#  my @assemble_execs = SelectFromDependants($working_dir, $test_build);

  my $diffs_file;
  for $diffs_file (keys %cvs_info)
  {

    # open gcov file
    if ( -e "$test_build/$diffs_file.gcov" )
    {
      # open output file
      open ( output_file, "> $working_dir/$diffs_file.counts.$tastyproject" )
	or die "Can't open $working_dir/$diffs_file.counts.$tastyproject";

      open( gcov_file, "< $test_build/$diffs_file.gcov" )
	or die "Can't open $test_build/$diffs_file.gcov";
    
      #  read line-by-line into an array
      my @gcov_array;
      while (<gcov_file>)
      {
	chomp;
	push @gcov_array, $_;
      }
      close(gcov_file);

      #maybe multiple sets of edits to a file
      my $i=0;
      while( $i <= (($#{ $cvs_info{$diffs_file} })-1) )
      {
	my $start_line = $cvs_info{$diffs_file}[$i];
	my $stop_line = $cvs_info{$diffs_file}[++$i];
	for (my $j=$start_line; $j<=$stop_line; $j++)
	{
	  print output_file "$gcov_array[$j]\n";
	}   
	++$i;
      }
    }   
  }
}
      

sub CleanGCOV
{
  my ($working_dir, $test_build) = @_;
  my @files_to_clean = ReadGCOVBuiltFile($working_dir);

  my $fil;
  my $res;
  foreach $fil (@files_to_clean)
  {
    $fil .= ".o";
    if ( -e "$test_build/$fil" )
    {
      my $cmmd_str = "rm $test_build/$fil";
      $res = system( $cmmd_str );

      if ( $res != 0 )
      {
	PrintMessage("\n$cmmd_str terminated with return code $res \n");
      }
    }
  }

  return $res == 0;
}

sub SelectFromDependants
{
  my @build_list;

  my ($working_dir, $build_dir) = @_;  
  my @dependants_list = ReadDependantsFile($working_dir);
  my @assemble_execs = ReadMakefile($build_dir);
  
  # if item in dependants_list is in assemble_execs, keep it
  my $deps;
  my $assm;
  foreach $deps (@dependants_list)
  {
    my $found_it = 0;
    foreach $assm (@assemble_execs)
    {
      my ($t,$h) = split(/\t/, $assm);
      $h =~ s/ //;
      my $dep_name_only = RemovePath($deps);

      if ($h eq $dep_name_only)
      {
	push @build_list, $deps;
	$found_it = 1;
      }
      last if ($found_it);
    }
  }

  return @build_list;
}

sub ReadMakefile
{
  my $test_build = $_[0];
 
  my $can_add = 0;
  my @assembly_execs;

  my $execs;
  open(X, "$test_build/Makefile")
    or die "Can't open $test_build/Makefile\n";
  while(<X>)
  {
    chomp;
    
    if ($can_add)
    {
      if ($_ =~ m/(\s*)(\w)(\s)\\/)
      {
	my ($exe, $postfx) = split(/\\/, $_);
	push @assembly_execs, $exe;
      }
      else
      {
	$can_add = 0;
      }
    }
  
    if ( $_ =~ m/ASSEMBLY_EXECS|UTILITY_EXECS|FINISHING_EXECS/ )
    {
      $can_add = 1;
    }
  }
  close(X);

  return @assembly_execs;
}

sub ReadDependantsFile
{

  my $working_dir = $_[0];
  my @make_list;
  open( dependants_file, "<$working_dir/tmp/dependants.all" ) 
    or die "Can't open $working_dir/tmp/dependants.all";
  while( <dependants_file> )
  {
    chomp;
    push @make_list, $_;
  }
  close(dependants_file);

  return @make_list;
}


# return last bit in a path
sub RemovePath
{
  my $full_bit = $_[0];

  my $end_bit;
  my @bits;
  if ( $full_bit !~ m/\// )
  {  
    $end_bit = $full_bit;
  }
  else
  {
    @bits = split(/\//, $full_bit);
    $end_bit = $bits[$#bits];    
  }
  return $end_bit;
}


#  grab each line in the cvs diffs file that belongs to the requested file
#  and is not a comment
sub ParseCVSDiffFile
{
  my ($working_dir, $test_build) = @_;

  open( cvs_diffs_file, "<$working_dir/tmp/cvs.diffs" ) 
    or die "Can't open $working_dir/tmp/cvs.diffs";

  my %modified_lines;
  my $indx;
  my $filename;
  my $start_line;
  my $stop_line;
  while ( <cvs_diffs_file> )
  {
    chomp;

    #"Index: " line gives us the filename
    if ( $_ =~ m/^Index:/)
    {
      ($indx, $filename) = split(/^Index: /, $_);
      $filename = RemovePath($filename);
    }

    elsif ( $_ =~ m/^\d+(,\d+)?[ac](\d+)(,\d+)?/ )
    {
      my ($s1, $s2) = split(/[ac]/, $_); 

      if ( $s2 =~ m/,/ )
      {
	($start_line, $stop_line) = split(/,/,$s2);
      }
      else
      {
	$start_line = $s2;
	$stop_line = $start_line;
      }

      push @{ $modified_lines{$filename} }, $start_line-1, $stop_line-1;

    }
  }
  close( cvs_diffs_file );
  return %modified_lines;
}

 

