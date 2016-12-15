#!/usr/bin/perl -w
package wGetOptions;

use Getopt::Long;
use Pod::Usage;
use base qw/Exporter/;

@EXPORT_OK = qw(wpod2usage wGetOptions);

$podstr = "";

sub wGetOptions {
    my @arg = @_;
    my %htypes = ( "s" => "STR", "i" => "INT", "f" => "FLT" );
    my $man = 0;
    my $help = 0;
    my @opts = ("help|?",\$help,"man",\$man);
    my @keys = qw(help man);
    my @types = ("","");
    my @shorts = ("Print out brief help message","Print the full documentation in man page style");
    my @longs = ("Print a help message and exits","Prints a manual page and exits upon typing 'q'");
    my @defaults = (0,0);
    my @isects = (0);  ## section indices
    my @tsects = ("General Options");  ## section titles
    my $main = "";
    for(my $i=0; $i < @_; ++$i) {
	if ( $arg[$i] =~ /^-/ ) {
	    if ( $arg[$i] =~ /^--/ ) {
		$arg[$i] =~ s/^-*//;
		push(@isects,$#keys+1);
		push(@tsects,$arg[$i]);
	    }
	    else {
		$arg[$i] =~ s/^-*//;
		$main = $arg[$i];
	    }
	}
	else {
	    my $opt = $arg[$i];
	    my ($ref,$short,$long) = @{$arg[$i+1]};
	    my ($key,$type) = split(/=/,$opt);
	    $long = $short unless ( defined($long) );
	    my $typestr = defined($type) ? $htypes{$type} : "";
	    my $default;
	    if ( ref($ref) eq "ARRAY" ) {
		$short = "(Multiples) $short";
		$long = "(Allows multiple values) $long";
		$default = join(" ",@{$ref});
	    }
	    else {
		if ( defined($type) ) {
		    $default = ${$ref};
		}
		else {
		    $default = ${$ref} ? "ON" : "OFF";
		}
	    }
	    push(@keys,$key);
	    push(@opts,$opt);
	    push(@opts,$ref);
	    push(@types,$typestr);
	    push(@shorts,$short);
	    push(@longs,$long);
	    push(@defaults,$default);
	    ++$i;
	}
    }

    my $ret = GetOptions(@opts);

    $podstr = "=pod\n\n=head1 NAME\n\n$0 - $main \n\n=head1 SYNOPSIS\n\n$0 [options]\n\n";
# General Options:\n";
#  -help             Print out brief help message [OFF]\n  -man              Print the full documentation in man page style [OFF]\n";
    my @values = ();
    for(my ($i,$j) = (0,0); $i < @keys; ++$i) {
	if ( ( $j <= $#isects ) && ( $i == $isects[$j] ) ) {
	    $podstr .= "\n $tsects[$j]:\n";
	    ++$j;
	}
	my $value;
	if ( ref($opts[$i+$i+1]) eq "ARRAY" ) {
	    $value = join(" ",@{$opts[$i+$i+1]});
	}
	else {
	    if ( $types[$i] ) {
		$value = ${$opts[$i+$i+1]};
	    }
	    else {
		$value = ${$opts[$i+$i+1]} ? "ON" : "OFF";
	    }
	}
	push(@values,$value);
	$podstr .= sprintf("  -%-17s%s [%s]\n","$keys[$i] $types[$i]",$shorts[$i],$values[$i]);
    }
    $podstr .= "\n=head1 OPTIONS\n\n=over 8\n\n=item B<-help>\n\nPrint a brief help message and exits\n\n=item B<-man>\n\nPrints the manual page and exits\n\n";
    for(my $i=0; $i < @keys; ++$i) {
	$podstr .= sprintf("=item -%s [%s]\n\n$longs[$i]\n\n","B<-$keys[$i] $types[$i]>",$values[$i]);
    }

    wpod2usage(-verbose => 1, -exitval => 1) if ( $help );
    wpod2usage(-verbose => 2) if ( $man );
    
    return $ret;
}

sub wpod2usage {
    ## Parse the input argument, same to what pod2usage does
    local($_) = shift;
    my %opts = ();
    ## Collect arguments
    if (@_ > 0) {
        ## Too many arguments - assume that this is a hash and
        ## the user forgot to pass a reference to it.
        %opts = ($_, @_);
    }
    elsif (!defined $_) {
	$_ = '';
    }
    elsif (ref $_) {
        ## User passed a ref to a hash
        %opts = %{$_}  if (ref($_) eq 'HASH');
    }
    elsif (/^[-+]?\d+$/) {
        ## User passed in the exit value to use
        $opts{'-exitval'} =  $_;
    }
    else {
        ## User passed in a message to print before issuing usage.
        $_  and  $opts{'-message'} = $_;
    }

    ## Temporarily write a POD document
    my $pid = $$;
    my @binpaths = split(/\//,$0);
    my $binname = pop(@binpaths);
    mkdir("/tmp/$pid");
    open(OUT,">/tmp/$pid/$binname") || die "Cannot open file\n";
    print OUT $podstr;
    close OUT;

    ## Modify the options to include -input, without exiting
    my $exitval = $opts{'-exitval'};
    $exitval = 0 unless ( defined($exitval) );
    $opts{'-exitval'} = "noexit";
    $opts{'-verbose'} = 0 unless defined($opts{'-verbose'});
    $opts{'-input'} = "/tmp/$pid/$binname";
 
    ## Call pod2usage function
    pod2usage(\%opts);
    #pod2usage({-exitval => "noexit", -input => "/tmp/$pid/$binname", -verbose => 2});
    unlink("/tmp/$pid/$binname");
    rmdir("/tmp/$pid");
    exit($exitval)  unless ($exitval eq 'noexit');
}
1;
