#!/usr/bin/perl

## Macd location (get from my github)
use lib '/home/oliver/GIT_REPOS/macd/lib';
use Macd;
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Find;
use strict;
use Data::Dumper;

## get these from environment variables

my $MODE=shift;


my $HOME = $ENV{HOME} || die "NO \$HOME variable set\n";
my $GRPATH = $ENV{GRPATH} || die "NO \$GRPATH variable set\n";
my $R_LIBS = $ENV{R_LIBS} || die "NO \$R_LIBS variable\n";
my $CHIGP_PATH = "$ENV{GRPATH}/CHIGP/";


## site specific conf file for macd - see http://github.com/ollyburren/macd
## for more details
my $grid_cnf = "$GRPATH/macd/example/ini/example.cnf";

my $r_lib_dir = $R_LIBS;
	
my $DRIVER = Macd::GRIDDriverFactory->instantiate('SGE',inifile => $grid_cnf);



###############
#CONFIGURATION#
###############
## where to mail to once finished.

my $BASE_DIR = "$GRPATH/astle/DATA/";

my $INDIR = "${BASE_DIR}gwas/";
my $OUTDIR = "${BASE_DIR}out/pmi/";
my $BLOCKFILE = "${BASE_DIR}support/0.1cM_regions.b37.bed";
my $W = 0.15;
my $DOTEST = 0; #IF set to true allows us to test the script by running only a few jobs

my $RSCRIPT="/home/oliver/bin/Rscript $GRPATH/astle/R/computePPI_nosdY.R";
my $LOG_DIR = "${BASE_DIR}/log/ppi/";

my $step = Macd::Step->new(
	logdir=> $LOG_DIR,
	driver=>$DRIVER,
	env_settings=>"R_LIBS=$r_lib_dir,GRPATH=$GRPATH"
);

if(! -e $OUTDIR){
	print "Making path $OUTDIR\n";
	make_path($OUTDIR);
}

my @INFILES;

find(sub{
		if(/common.tsv$/){
			(my $stub = $_)=~s/\.tsv//;
			push @INFILES,$stub;
		}
},($INDIR));

#die Dumper(\@INFILES);

my $addcount;
foreach my $i(@INFILES){
	my $ofile = "$OUTDIR/${i}.pmi";
	next if -e $ofile;
	my @param = "$RSCRIPT";
	push @param, "--in.file=$INDIR/$i.tsv";
	push @param, "--out.file=$ofile";
	push @param, "--W=$W";
	push @param, "--block.file=$BLOCKFILE";
	my $cmd = join(" ",@param);
	#die "$cmd";
	my $job =  Macd::Step::Job->new(command=>$cmd);
	$step->add_job($job);
	$addcount++;
	last if $DOTEST;
}

print "Added $addcount jobs\n";

if($step->execute()){
	print "Step submitted successfully\n";
	## we can hold up prog execution as follows
	## in case we have a step below that requires the output
	$step->wait_on_complete();
	print "Step completed successfully\n";
}else{
	print "Error submitting step\n";
}
