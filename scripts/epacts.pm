package epacts;

use strict;
use base qw/Exporter/;
use lib "$FindBin::Bin";
use File::Basename;

## Variables and methods shared across the package
our @EXPORT = qw(%hszchrs @chrs @szchrs @cumszchrsMb %ichrs);
our @EXPORT_OK = qw(parsePheno getMosixCmd schr2nchr vcfSampleIDs vcfSampleIndex readPedVcf readPedVcfMulti readPedKinLabel $binR $binRscript $binrm $binmake $binzcat $bincat $binhead $binmv $bincut $bingrep $binawk $binpfbtops $bingnuplot $binepstopdf $binsort $defaultfasta $xLabel $yLabel $xStart $xStop installPackages tofpos fromfpos forkExecWait initRef intervalsByUnit intervalsByBED readBED);

our $epactsdir = dirname($FindBin::Bin);
our $datadir = "$epactsdir/share/EPACTS";
our $defaultfasta = "$datadir/human_g1k_v37.fasta";

our $binR = "R";
our $binRscript = "Rscript"; #$binRscript = "R CMD BATCH --slave --no-save --no-restore";
our $binrm = "rm";
our $binmake = "make";
our $binzcat = "zcat";
our $bincat = "cat";
our $binhead = "head";
our $binmv = "mv";
our $bincut = "cut";
our $bingrep = "grep";
our $binawk = "awk";
our $binsort = "sort";
our $binpfbtops = "pfbtops";
our $bingnuplot = "gnuplot";
our $binepstopdf = "$epactsdir/bin/epstopdf";
our $xLabel = "X";
our $yLabel = "Y";
our $mtLabel = "MT";
our $xStart = 2699520;
our $xStop = 154931044;

our %hszchrs = ();
our @chrs   = ();
our @szchrs = ();
our @cumszchrsMb = ();
our %ichrs   = ();
#my  %hchrs   = ();

#BEGIN {
## Variables below are hard-coded based on GRCh37
#    @chrs = (1..22,"X","Y","MT");
#    @szchrs = qw(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560 59373566 16569);
#    %ichrs = ();
#    
#    @cumszchrsMb = (0);
#    for(my $i=0; $i < @chrs; ++$i) {
#	push(@cumszchrsMb,$szchrs[$i]/1e6+$cumszchrsMb[$i]);
#	$ichrs{$chrs[$i]} = $i;
#    }
#}

## parse phenotype value, resolving missing values
sub parsePheno {
    my ($p,$missing) = @_;
    if ( $p eq $missing ) {
	return "NA";
    }
    elsif ( $p =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
	return $p;
    }
    else {
	return "NA";  ## non-numeric values are converted to NA
    }
}

## convert command to mosix command
sub getMosixCmd {
    my ($cmd,$nodes) = @_;
    return "mosbatch -E/tmp -i -j$nodes sh -c '$cmd'";
}

## convert string chromosome to integer (compatible to PLINK)
sub schr2nchr {
    my $chr = shift;
    $chr =~ s/^chr// if ( substr($chr,0,3) eq "chr" );
    if ( $chr =~ /^\d+$/ ) { return $chr; }
    elsif ( $chr eq "X" ) { return 23; }
    elsif ( $chr eq "Y" ) { return 24; }
    elsif ( $chr eq "MT" ) { return 25; }
}

sub vcfSampleIDs {
    my $vcf = shift;
    my @F = ();
    if ( -s "$vcf.tbi" || -s "$vcf.csi" || -s "$vcf.s1r" ) {
	@F = split(/[\t\r\n]+/,`$epactsdir/bin/parse-sample-ids $vcf`); #tabix -H $vcf | tail -1`);
    }
    else {
	@F = split(/[\t\r\n]+/,`$binzcat $vcf | $binhead -1000 | $bingrep ^#CHROM`);
        splice(@F,0,9);
    }
    return @F;
}

sub vcfSampleIndex {
    my ($vcf,$indf) = @_;
    my @ids = &vcfSampleIDs($vcf);
    my @idxs = ();
    if ( $indf eq "" ) {
	for(my $i=0; $i < @ids; ++$i) {
	    push(@idxs,$i);
	}
    }
    else {
	my %hids = ();
	if ( $indf ) {
	    open(IN,$indf) || die "Cannot open file $indf\n";
	    while(<IN>) {
		my ($id) = split;
		$hids{$id} = 1;
	    }
	    close IN;
	}
	for(my $i=0; $i < @ids; ++$i) {
	    push(@idxs,$i) if ( defined($ids[$i]) );
	}
    }
    return (@idxs);
}

sub readPedVcf {
    my ($epactsdir,$ped,$vcf,$missing,$pheno,$rcovs,$rcondsnps,$field,$pass) = @_;
    $rcovs = [] unless ( defined($rcovs) );
    $rcondsnps = [] unless ( defined($rcondsnps) );
    $missing = "NA" unless ( defined($missing) );
    $field = "GT" unless ( $field );
    $pass = "" unless ( ( defined($pass) ) && $pass );

    my %hVcfIds = ();
    my @covs = @{$rcovs};
    my @condsnps = @{$rcondsnps};

    my @vcfIds = split(/[\s\t\r\n]+/,`$epactsdir/bin/parse-sample-ids $vcf`); #tabix -H $vcf | tail -1`);
    #splice(@vcfIds,0,9);
    for(my $i=0; $i < @vcfIds; ++$i) { $hVcfIds{$vcfIds[$i]} = $i; }
    
    my @datIds = ();
    my $iphe;
    my @icovs = ();
    
    open(PED,$ped) || die "Cannot open PED file $ped\n";
    my $dat = $ped; 
    $dat =~ s/\.ped$/\.dat/;
    if ( -s $dat ) { ## check if dat file exist
	open(DAT,$dat) || die "Cannot open DAT file $dat\n";
	@datIds = qw(FAM_ID IND_ID FAT_ID MOT_ID SEX);
	for(my $i=$#datIds+1;<DAT>;++$i) {
	    my ($symbol,$name) = split(/[\s\t\r\n]+/);
	    push(@datIds,$name);
	}
	close DAT;
    }
    else {
	my @F = split(/[\s\r\t\n]/,`$binhead -1 $ped`);
	if ( $F[0] =~ /^#/ ) {
	    $F[0] =~ s/^#//;
	    @datIds = @F;
	}
    }
    
    if ( $#datIds < 0 ) {  ## NO HEADER available
	$iphe = 5;
	#if ( $allcov ) {
	#    my @F = split(/[\s\t\r\n]+/,`$binhead -1 $ped`);
	#    for(my $i=6; $i < @F; ++$i) {
	#	push(@icovs,$i);
	#    }
	#}
	if ( $#covs >= 0 ) {
	    die "FATAL ERROR: Cannot find $covs[0] in $ped";	    
	}
    }
    else {  ## HEADER is available
	my %hcols = ();
	for(my $i=0; $i < @datIds; ++$i) {
	    $hcols{$datIds[$i]} = $i;
	}
	
	if ( $pheno ) {
	    $iphe = $hcols{$pheno};
	    die "Cannot find $pheno in $ped" unless defined($iphe);
	}
	else {
	    $iphe = 5;
	}
	
	#if ( $allcov ) {
	#    for(my $i=5; $i < @datIds; ++$i) {
	#	push(@icovs,$i) unless ( $i == $iphe );
	#    }
	#}
	#else {
	    for(my $i=0; $i < @covs; ++$i) {
		my $icov = $hcols{$covs[$i]};
		die "FATAL ERROR: Cannot find $covs[$i] in $ped" unless defined($iphe);
		push(@icovs,$icov);
	    }
	#}
    }

    ## Parse VCF and put condition SNP info
    my @vcfCondGenos = ();
    foreach my $condsnp (@condsnps) {
	my @genos = &vcfExtractMarkerGenotypes($vcf,$condsnp,$field,$pass);
	if ( $#genos < 0 ) {
	    print STDERR "WARNING: Cannot find $condsnp from VCF file... Skipping..\n";
	    next;
	}

	die "# of VCF IDs and genotypes do not match $#vcfIds vs $#genos\n" unless ( $#vcfIds == $#genos );
	my $ac = 0;
	my $an = 0;
	for(my $i=0; $i < @vcfIds; ++$i) { 
	    if ( $genos[$i] ne "NA" ) {
		$ac += $genos[$i];
		$an += 2;
	    }
	}

	## Apply mean imputation for missing genotypes
	for(my $i=0; $i < @genos; ++$i) {
	    $genos[$i] = (( $an == 0 ) ? 0 : sprintf("%.4lf",($ac/$an))) if ( $genos[$i] eq "NA" );
	}
	push(@vcfCondGenos,\@genos);
    }
    
    my %hPhes = ();
    my %hCovs = ();
    my $nPhes = 0;
    for(my $i=0;<PED>;++$i) {
	next if ( /^#/ );
	my @F = split(/[\s\t\r\n]+/);
	my $id = $F[1];
	
	die "ERROR: # of columns in $ped (".($#F+1).") do not match with previous lines (".($#datIds+1).") at line $. Note that any whitespace will be considered as delimiters\n" unless ($#F == $#datIds);
	
	if ( defined($hVcfIds{$id}) ) {
	    my $p = defined($pheno) ? &parsePheno($F[$iphe],$missing) : "1";
	    next if ( $p eq "NA" );
	    $hPhes{$id} = $p;
	    my @c = ();
	    for(my $j=0; $j < @icovs; ++$j) {
		push(@c,&parsePheno($F[$icovs[$j]],$missing));
		#die "ERROR: Missing covariate value is detected. Currently EPACTS won't run correctly with missing covariates (bug fix TBA)\n" if ( $c[$#c] eq $missing );
	    }
	    my $ivcf = $hVcfIds{$id};
	    for(my $j=0; $j < @vcfCondGenos; ++$j) {
		push(@c,$vcfCondGenos[$j]->[$ivcf]);
	    }
	    $hCovs{$id} = \@c;
	    ++$nPhes;
	}
    }
    
    die "ERROR: No overlapping IDs between VCF and PED file. Cannot proceed.\n" if ( $nPhes == 0 );
    
## check if binary phenotypes and convert if needed
    my %valPhes = ();
    foreach my $id (keys %hPhes) {
	my $val = $hPhes{$id};
	next if ( $val eq $missing );
	$valPhes{$val} = 0 unless ( defined($valPhes{$val}) );
	++$valPhes{$val};
    }
    my @uniqValPhes = sort {$a <=> $b} keys %valPhes;
    my $isBinary = 0;
    if ( $#uniqValPhes > 1 ) {
	print STDERR "Detected phenotypes with more than 3 unique values -- considering as quantitative phenotypes.\n";
    }
    elsif ( $#uniqValPhes == 1 ) {
	print STDERR "Detected phenotypes with 2 unique values - $uniqValPhes[0] and $uniqValPhes[1] - considering them as binary phenotypes... re-encoding them into 1 and 2\n";
	foreach my $id (keys %hPhes) {
	    my $val = $hPhes{$id};
	    unless ( $val eq $missing ) {
		if ( $val == $uniqValPhes[0] ) {
		    $hPhes{$id} = 1;
		}
		elsif ( $val == $uniqValPhes[1] ) {
		    $hPhes{$id} = 2;
		}
		else {
		    die "Cannot recognize binary phenotype value $val\n";
		}
	    }
	}	    
	$isBinary = 1;
    }
    elsif ( defined($pheno) ) {
	die "ERROR: Phenotypes has only one or less unique values (@uniqValPhes)";
    }

    return(\@vcfIds,\%hPhes,\%hCovs,$isBinary);
}

sub readPedKinLabel {
    my ($epactsdir,$ped,$evec,$label) = @_;

    my %hVcfIds = ();

    my @vcfIds = split(/[\s\t\r\n]+/,`cut -f 1 $evec | grep -v ^#`);
    for(my $i=0; $i < @vcfIds; ++$i) { $hVcfIds{$vcfIds[$i]} = $i; }

    my @datIds = ();
    my $iphe;
    
    open(PED,$ped) || die "Cannot open PED file $ped\n";
    my $dat = $ped; 
    $dat =~ s/\.ped$/\.dat/;
    if ( -s $dat ) { ## check if dat file exist
	open(DAT,$dat) || die "Cannot open DAT file $dat\n";
	@datIds = qw(FAM_ID IND_ID FAT_ID MOT_ID SEX);
	for(my $i=$#datIds+1;<DAT>;++$i) {
	    my ($symbol,$name) = split(/[\s\t\r\n]+/);
	    push(@datIds,$name);
	}
	close DAT;
    }
    else {
	my @F = split(/[\s\r\t\n]/,`$binhead -1 $ped`);
	if ( $F[0] =~ /^#/ ) {
	    $F[0] =~ s/^#//;
	    @datIds = @F;
	}
    }
    
    if ( $#datIds < 0 ) {  ## NO HEADER available
	$iphe = 5;
    }
    else {  ## HEADER is available
	my %hcols = ();
	for(my $i=0; $i < @datIds; ++$i) {
	    $hcols{$datIds[$i]} = $i;
	}
	
	if ( $label ) {
	    $iphe = $hcols{$label};
	    die "Cannot find $label in $ped" unless defined($iphe);
	}
	else {
	    $iphe = 5;
	}
    }

    my %hPhes = ();
    #my $nPhes = 0;
    for(my $i=0;<PED>;++$i) {
	next if ( /^#/ );
	my @F = split(/[\s\t\r\n]+/);
	my $id = $F[1];
	
	die "ERROR: # of columns in $ped (".($#F+1).") do not match with previous lines (".($#datIds+1).") at line $. Note that any whitespace will be considered as delimiters\n" unless ($#F == $#datIds);
	
	if ( defined($hVcfIds{$id}) ) {
	    my $p = $F[$iphe];
	    $hPhes{$id} = $p;
	}
    }
    
    #die "ERROR: No overlapping IDs between VCF and PED file. Cannot proceed.\n" if ( $nPhes == 0 );
    
    return(\@vcfIds,\%hPhes);
}

sub tofpos {
    my ($chr,$pos) = @_;
    return sprintf("%d.%09d",$ichrs{$chr},$pos);
}

sub fromfpos {
    my ($chr,$pos) = split(/\./,$_[0]);
    $chr = $chrs[$chr];
    $pos =~ s/^0*//;
    return ($chr,$pos);
}

sub forkExecWait {
    my $cmd = shift;
    print "forkExecWait(): $cmd\n";
    my $kidpid;
    if ( !defined($kidpid = fork()) ) {
	die "Cannot fork: $!";
    }
    elsif ( $kidpid == 0 ) {
	exec($cmd);
	die "Cannot exec $cmd: $!";
    }
    else {
	waitpid($kidpid,0);
    }
    die "Error in running $cmd. Exit code is ".($? >> 8)."\n" unless ( $? >> 8 == 0 );
}

sub readPedVcfMulti {
    my ($epactsdir,$ped,$vcf,$missing,$rphes,$rcovs,$rcondsnps,$field) = @_;
    $rphes = [] unless ( defined($rphes) );
    $rcovs = [] unless ( defined($rcovs) );
    $rcondsnps = [] unless ( defined($rcondsnps) );
    $missing = "NA" unless ( defined($missing) );

    my %hVcfIds = ();
    my @phes = @{$rphes};
    my @pnames = ();
    my @covs = @{$rcovs};
    my @condsnps = @{$rcondsnps};

    #print STDERR "foo $epactsdir\n";

    my @vcfIds = split(/[\s\t\r\n]+/,`$epactsdir/bin/parse-sample-ids $vcf`); #tabix -H $vcf | tail -1`);
    #splice(@vcfIds,0,9);
    for(my $i=0; $i < @vcfIds; ++$i) { $hVcfIds{$vcfIds[$i]} = $i; }

    #print STDERR "bar\n";
    
    ## read header information from either first line of PED of separate .dat file
    my @datIds = ();
    my @iphes = ();
    my @icovs = ();
    
    open(PED,$ped) || die "Cannot open PED file $ped\n";
    my $dat = $ped; 
    $dat =~ s/\.ped$/\.dat/;
    if ( -s $dat ) { ## check if dat file exist
	open(DAT,$dat) || die "Cannot open DAT file $dat\n";
	@datIds = qw(FAM_ID IND_ID FAT_ID MOT_ID SEX);
	for(my $i=$#datIds+1;<DAT>;++$i) {
	    my ($symbol,$name) = split(/[\s\t\r\n]+/);
	    push(@datIds,$name);
	}
	close DAT;
    }
    else {
	my @F = split(/[\s\r\t\n]/,`$binhead -1 $ped`);
	if ( $F[0] =~ /^#/ ) {
	    $F[0] =~ s/^#//;
	    @datIds = @F;
	}
    }
    
    die "Cannot find header line of PED or separate .dat file\n" if ( $#datIds < 0 );
    
    ## make a hash between trait name and column index
    my %hcols = ();
    for(my $i=0; $i < @datIds; ++$i) {
	$hcols{$datIds[$i]} = $i;
	#$hicols{$i} = 1;
    }

    my %hisels = ();
    ## determine the column indices of covariates first
    for(my $i=0; $i < @covs; ++$i) {
	my $icov = $hcols{$covs[$i]};
	die "FATAL ERROR: Cannot find $covs[$i] in $ped" unless defined($icov);
	push(@icovs,$icov);
	$hisels{$icov} = 1;
    }

    ## determine the column indices of the phenotypes
    if ( $#phes < 0 ) { # no phenotype selected - select the rest of phenotypes
	for(my $i=5; $i < @datIds; ++$i) {
	    unless ( defined($hisels{$i}) ) {
		push(@iphes,$i);
		push(@phes,$datIds[$i]);
	    }
	}
    }
    else {
	for(my $i=0; $i < @phes; ++$i) {
	    my $iphe = $hcols{$phes[$i]};
	    die "FATAL ERROR: Cannot find $phes[$i] in $ped" unless defined($iphe);

	    if ( defined($hisels{$iphe}) ) {
		die "FATAL ERROR: phenotype $phes[$i] is also included in the covariates\n";
	    }
	    else {
		push(@iphes,$iphe);
	    }
	}
    }

    ## Parse VCF and put condition SNP info
    my @vcfCondGenos = ();
    foreach my $condsnp (@condsnps) {
	my @genos = &vcfExtractMarkerGenotypes($vcf,$condsnp,$field);
	if ( $#genos < 0 ) {
	    print STDERR "WARNING: Cannot find $condsnp from VCF file... Skipping..\n";
	    next;
	}

	die "# of VCF IDs and genotypes do not match\n" unless ( $#vcfIds == $#genos );
	my $ac = 0;
	my $an = 0;
	for(my $i=0; $i < @vcfIds; ++$i) { 
	    $ac += $genos[$i];
	    $an += 2;
	}

	## Apply mean imputation for missing genotypes
	for(my $i=0; $i < @genos; ++$i) {
	    $genos[$i] = (( $an == 0 ) ? 0 : sprintf("%.4lf",($ac/$an))) if ( $genos[$i] eq "NA" );
	}
	push(@vcfCondGenos,\@genos);
    }
    
    ## Extract the phenotype and covariate data from PED file
    my %hPhes = ();
    my %hCovs = ();
    my $nPhes = 0;
    for(my $i=0;<PED>;++$i) {
	next if ( /^#/ );
	my @F = split(/[\s\t\r\n]+/);
	my $id = $F[1];
	
	die "ERROR: # of columns in $ped (".($#F+1).") do not match with previous lines (".($#datIds+1).") at line $. Note that any whitespace will be considered as delimiters\n" unless ($#F == $#datIds);
	
	if ( defined($hVcfIds{$id}) ) {
	    my @p = ();
	    for(my $j=0; $j < @iphes; ++$j) {
		push(@p,&parsePheno($F[$iphes[$j]],$missing));
		#die "ERROR: Missing phenotype value is detected in individual $id, phenotype $rphes->[$j]. Currently EPACTS won't run with missing phenotypes\n" if ( $p[$#p] eq $missing );
	    }
	    $hPhes{$id} = \@p;
	    
	    my @c = ();
	    for(my $j=0; $j < @icovs; ++$j) {
		push(@c,&parsePheno($F[$icovs[$j]],$missing));
		#die "ERROR: Missing covariate value is detected in individual $id, covariate $rcovs->[$j]. Currently EPACTS won't run correctly with missing covariates\n" if ( $c[$#c] eq $missing );
	    }
	    my $ivcf = $hVcfIds{$id};
	    for(my $j=0; $j < @vcfCondGenos; ++$j) {
		push(@c,$vcfCondGenos[$j]->[$ivcf]);
	    }
	    $hCovs{$id} = \@c;
	    ++$nPhes;
	}
    }
    
    die "ERROR: No overlapping IDs between VCF and PED file. Cannot proceed.\n" if ( $nPhes == 0 );
    return(\@vcfIds,\%hPhes,\%hCovs,\@phes);
}

sub vcfExtractMarkerGenotypes {
    my ($vcf,$markerId,$field,$pass) = @_;
    my @L = ();
    open(IN,"$epactsdir/bin/vcfast convert --vcf $vcf --marker-id $markerId --field $field ".($pass ? "" : "--ignoreFilter ")."--out - | grep -v ^##|") || die "Cannot execute vcfast\n";
    while(<IN>) {
	#print $_;
	my @F = split(/[\t\r\n]/);
	push(@L,\@F);
    }
    close IN;
    die "ERROR in accessing $vcf at marker $markerId\n" if ($? >> 8);
    if ( $#L < 0 ) {
	die "ERROR in reading $vcf at marker $markerId. No headers are found\n";
    }
    elsif ( $#L == 0 ) {
	my @empty = ();
	return @empty;
    }
    else {
	for(my $i=1; $i < @L; ++$i) {
	    if ( ( defined($L[$i]->[0]) ) && ( substr($L[$i]->[0],0,length($markerId)) eq $markerId ) ) {
		shift(@{$L[$i]});
		return @{$L[$i]};
	    }
	}
	my @empty = ();
	return @empty;
    }
}

sub initRef {
    my $ref = shift;

    unless ( %hszchrs ) {
	$ref = $defaultfasta unless ( ( defined($ref) ) && ( $ref ) );
	open(IN,"$ref.fai") || die "Cannot open file $ref.fai\n";
	my $cumbase = 0;

	my @autoChrs = ();
	my $chrPrefix = -1;
	## %hszchrs contains mapping between all chromosomes and the line
	for(my $nchr=0; <IN>; ++$nchr) {
	    my ($chrom,$base,$startbyte,$basesperline,$bytesperline) = split;
	    $hszchrs{$chrom} = [$startbyte,$basesperline,$bytesperline,$base,$nchr,$cumbase,$cumbase+$base];
	    $cumbase += $base;	    	    	    
	    if ( $chrom =~ /^chr(\d+)$/ ) {
		push(@autoChrs,$1);
		if ( $chrPrefix == 0 ) {
		    die "ERROR: Chromosomes with chr prefix and without chr prefix observed in $ref.fai";
		}
		elsif ( $chrPrefix == -1 ) {
		    $chrPrefix = 1;
		}
	    }
	    elsif ( $chrom =~ /^(\d+)$/ ) {
		push(@autoChrs,$1);
		if ( $chrPrefix == 1 ) {
		    die "ERROR: Chromosomes with chr prefix and without chr prefix observed in $ref.fai";
		}
		elsif ( $chrPrefix == -1 ) {
		    $chrPrefix = 0;
		}		
	    }
	    elsif ( ( $chrom eq "X" ) || ( $chrom eq "chrX" ) ) {
		## Hard-coded start, stop regions
		$xLabel = $chrom;
		if ( $base == 155270560 ) {
		    $xStart = 2699520;
		    $xStop  = 154931044;
		}
		elsif ( $base == 156040895 ) {
		    $xStart = 2781479;
		    $xStop  = 15570138;
		}
		else {
		    die "ERROR: The chromosome X size $base is not recognized by neither GRCh37 (hg19) nor GRCh38 (hg38). If you are using other genome builds, modify initRef function in epacts.pm to change PAR region info manually\n";
		}
	    }
	    elsif ( ( $chrom eq "Y" ) || ( $chrom eq "chrY") ) {
		$yLabel = $chrom;
	    }
	    elsif ( ( $chrom eq "MT" ) || ( $chrom eq "chrM" ) || ( $chrom eq "chrMT" ) ) {
		$mtLabel = $chrom;
	    }
	}

	die "Error in parsing FASTA index file $ref.fai\n" if ( ( $chrPrefix < 0 ) || ( ! $xLabel ) || ( ! $yLabel ) || ( ! $mtLabel ) );

	@chrs = ();
	@szchrs = ();
	@cumszchrsMb = (0);
	my @sortedAutoChrs = sort {$a <=> $b} @autoChrs;
	foreach my $c (@sortedAutoChrs,"X","Y") {
	    push(@chrs, ($chrPrefix == 1) ? "chr$c" : $c);
	    push(@szchrs, $hszchrs{$chrs[$#chrs]}->[3]);
	    my $newcumszMb = $szchrs[$#szchrs]/1e6 + $cumszchrsMb[$#cumszchrsMb];
	    push(@cumszchrsMb, $newcumszMb);
	    $ichrs{$chrs[$#chrs]} = $#chrs;
	}
	close IN;
	#@open(FASTA,$ref) || die "Cannot open file $ref\n";
    }
}

sub zopen {
    my $fn = shift;
    my $reg = shift;
    my $fh = FileHandle->new;
    if ( $fn =~ /\.gz$/ ) {
	die "Cannot open file $fn\n" unless ( -s $fn );
	if ( defined($reg) && ( $reg ) ) {
	    die "Cannot parse $reg\n" unless ( $reg =~ /^\S+:\d+(-\d+)?/ );
	    die "Cannot open file $fn.tbi\n" unless ( -s "$fn.tbi" );
	    open($fh,"$epactsdir/bin/tabix -h $fn $reg |");
	}
	else {
	    open($fh,"$binzcat $fn|");
	}
    }
    else {
	die "Cannot parse region $reg in a plain text\n" if ( defined($reg) && ( $reg ) );
	if ( $fn eq "-" ) {
	    $fh = *STDIN;
	}
	else {
	    open($fh,$fn) || die "Cannot open $fn\n";
	}
    }
    return $fh;
}

sub readBED {
    my $fn = shift;
    my $fh = zopen($fn);
    my %hbed = ();
    while(<$fh>) {
	my ($chr,$beg0,$end1) = split;
	die "Invalid BED interval $chr:$beg0-$end1 in $fn\n" if ( $beg0 >= $end1 );
	if ( defined($hbed{$chr}) ) {
	    my $r = $hbed{$chr};
	    my $nr = $#{$r};
	    die "Unsorted or overlapping BED interval $chr:".join("\t",@{$r->[$nr]})." and $chr:$beg0-$end1\n" if ( $r->[$nr]->[1] > $beg0 );
	    push(@{$hbed{$chr}},[$beg0,$end1]);
	}
	else {
	    $hbed{$chr} = [ [$beg0,$end1] ];
	}
    }
    return \%hbed;
}

sub intervalsByBED {
    my ($bedf,$rchrs,$rszchrs,$chrom,$bprange) = shift;
    my $rBED = &readBED($bedf);
    my @intervals = ();    
    for(my $i=0; $i <= $#{$rchrs}; ++$i) {  ## for each chromosome..
	my @chrIntervals = ();

	if ( defined($chrom) && ( $chrom ) && ( $rchrs->[$i] ne $chrom ) ) {
	    push(@intervals,\@chrIntervals);
	    next;
	}
	
	my $chr = $rchrs->[$i];
	my $r = $rBED->{$chr};
	if ( defined($r) ) {
	    foreach my $i (@{$r}) {
		my $beg = $i->[0]+1;
		my $end = $i->[1];
		
		## further subset the region if possible
		if ( $bprange ) {
		    my ($b,$e) = split(/-/,$bprange);
		    ## no overlap if
		    ## $b < $e < $beg < $end
		    ## $beg < $end < $b < $e
		    next if ( ( $e < $beg ) || ( $end < $b ) );
		    ## $beg < $e,  $b < $end is guaranteed
		    ## the interval is max($b,$beg) to min($e,$end)
		    $beg = $b if ( $b > $beg );
		    $end = $e if ( $e < $end );
		}
		
		push(@chrIntervals, [$beg,$end]);
	    }
	}
	push(@intervals,\@chrIntervals);
    }
    return \@intervals;
}

sub intervalsByUnit {
    my ($unit,$rchrs,$rszchrs,$rexbed,$chrom,$bprange) = @_;
    my @intervals = ();
    for(my $i=0; $i <= $#{$rchrs}; ++$i) {  ## for each chromosome..
	my @chrIntervals = ();
	
	if ( defined($chrom) && ( $chrom ) && ( $rchrs->[$i] ne $chrom ) ) {
	    push(@intervals,\@chrIntervals);
	    next;
	}
	
	my $chr = $rchrs->[$i];
	my $szchr = $rszchrs->[$i];

	my @exbeds = ();
	## avoid BED region if possible
	## starts from [$beg,$end]
	if ( defined($rexbed) && defined($rexbed->{$chr}) ) {
	    @exbeds = @{$rexbed->{$chr}}; 

	    for(my $j=0; $j < @exbeds; ++$j) {
		die "Invalid BED interval $chr:".join("-",@{$exbeds[$j]})."\n" unless ( $exbeds[$j]->[0] < $exbeds[$j]->[1] );
		die "Unsorted or overlapping BED interval $chr:".join("-",@{$exbeds[$j-1]})." and $chr:".join("-",@{$exbeds[$j]})."\n" if ( ( $j > 0 ) && ( $exbeds[$j-1]->[1] < $exbeds[$j]->[0] ) );
	    }
	}
	
	for(my $beg=1; $beg < $szchr; $beg += $unit) {
	    my $end = $beg + $unit - 1;
	    $end = $szchr if ( $end > $szchr );

	    ## $chr:$beg-$end

	    ## further subset the region if possible
	    if ( $bprange ) {
		my ($b,$e) = split(/-/,$bprange);
		## no overlap if
		## $b < $e < $beg < $end
		## $beg < $end < $b < $e
		next if ( ( $e < $beg ) || ( $end < $b ) );
		## $beg < $e,  $b < $end is guaranteed
		## the interval is max($b,$beg) to min($e,$end)
		$beg = $b if ( $b > $beg );
		$end = $e if ( $e < $end );
	    }

	    if ( $#exbeds >= 0 ) {
		## find intervals that overlaps
		my $beg0 = $beg-1;
		my $end1 = $end;
		my $j = 0;
		while ( ( $j <= $#exbeds ) && ( $exbeds[$j]->[1] <= $beg0 ) ) {
		    ++$j;
		}
		next if ( $j > $#exbeds );
		next if ( $exbeds[$j]->[0] >= $end1 );

		## there is overlap
		while( $beg0 < $end1 ) {
		    if ( $exbeds[$j]->[0] > $beg0 ) {
			push(@chrIntervals,[$beg0+1,$exbeds[$j]->[0]]);
		    }
		    $beg0 = $exbeds[$j]->[1];
		    ++$j;
		    last if ( ( $j > $#exbeds ) || ( $exbeds[$j]->[0] >= $end1 ) );
		}
		push(@chrIntervals,[$beg0+1,$end1]) if ( $beg0 < $end1 )
	    }
	    else {
		push(@chrIntervals,[$beg,$end]);
	    }
	}
	push(@intervals,\@chrIntervals);
    }
    return \@intervals;
}

1;
