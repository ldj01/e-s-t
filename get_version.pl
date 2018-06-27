#!/usr/bin/perl -w
#------------------------------------------------------------------------------
# Name: get_version.pl
#
# Description:	Generates a version string, incorporating
# relevant Git information into it.
#
#------------------------------------------------------------------------------
use strict;
use Cwd;
use File::Basename;

# Forward declarations
sub ShowUsage();
sub GitGetVersion();

# Option variables
my $help = '';

# Main variables
my $source_dir;          # main source code directory; home of version.txt
my $revision = '';
my $is_modified = ''; # Whether there are pending changes under the working
                      # path.
my $output = '';

# Verify that we have the required input path to the main build directory.
ShowUsage() if ($#ARGV < 0);
$source_dir = $ARGV[0];

# See if we're in a Git repo.  If not, then read the version info from
# version.txt.
my $git_dir = `git rev-parse --git-dir 2>> /dev/null`;
if ($? != 0)
{
    if (-f "$source_dir/version.txt")
    {
        chomp($output = `cat $source_dir/version.txt`);
    }
    else
    {
        print STDERR "Warning: $source_dir/version.txt doesn't exist.\n";
        $output = 'Unknown';
    }

    print "$output\n";
    exit(0);
}

# We're in a Git repository, so get the version using Git commands.
$output = GitGetVersion();

# Store the version in version.txt in the main source directory.
open(VERSION, ">$source_dir/version.txt")
    or die "Couldn't create $source_dir/version.txt: $!.";
print VERSION "$output\n";
close(VERSION) or die "Couldn't close $source_dir/version.txt.";

print "$output\n";
exit(0);


#------------------------------------------------------------------------------
# Name: GitGetVersion()
# Get the tag or branch name of currently checked out item.  
# Note that if the name cannot be determined, The name "Unknown" is returned.
#------------------------------------------------------------------------------
sub GitGetVersion()
{
    # Attempt to get branch name
    my $tname = `git rev-parse --abbrev-ref HEAD 2>> /dev/null`;
    return ("Unknown") if ($? != 0);

    # If the name is 'HEAD' then we have a tag instead, so return the
    # real tag name.
    chomp $tname;
    if ($tname eq "HEAD")
    {
        chomp($tname = `git describe --tags`);
        return $tname;
    }

    # Not a tag, so get the revision.
    chomp($revision = `git rev-parse --short HEAD 2>> /dev/null`);

    # No revision, so just return the branch name.
    return $tname if (not $revision);

    $tname .= "-$revision";

    # Find out if there are any local mods
    chomp($is_modified = `git status -s -uno 2>> /dev/null`);

    # If there are local mods, tack on "M"
    $tname .= "M" if ($is_modified ne '');

    return $tname;
}

#------------------------------------------------------------------------------
# Name: ShowUsage()
#
# Description: Shows the help or usage directives for the program.
#------------------------------------------------------------------------------
sub ShowUsage()
{
    print << "EOF";
Create a version for build, using Git if in repository.

  get_version.pl source_dir
      source_dir    path to the main source directory
EOF
    exit(1);
}
