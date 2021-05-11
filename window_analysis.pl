use strict;

use Getopt::Long qw( :config posix_default bundling no_ignore_case );;
use Scalar::Util qw(looks_like_number);

my $file = '';
my $chr = 1;
my $start = 2;
my $end = 3;
my $window = 1000;
my $start_ind = '';
my $end_ind = '';
my $control_min = 24;
my $control_max = 23;
my $control_q1 = 25;
my $control_q9 = 26;
my $control_mean = 23;
my $control_mean_thres = 0.2;
my $b_thres = 0.1;
my $n_thres = 2;
my $mode = "BOTH";
my $outfile = "STDOUT";


help() if ( @ARGV < 1);
GetOptions('file|f=s' => \$file,
           'chr|c=i' => \$chr,
           'start|s=i' => \$start,
           'end|e=i' => \$end,
           'window|w=i' => \$window,
           'start_ind|a=i' => \$start_ind,
           'end_ind|b=i' => \$end_ind,
           'c_min|n=i' => \$control_min,
           'c_max|x=i' => \$control_max,
           'c_q1|1=i' => \$control_q1,
           'c_q9|9=i' => \$control_q9,
           'c_mean|M=i' => \$control_mean,
           'c_mean_thres|u=f' => \$control_mean_thres,
           'beta_thres|t=f' => \$b_thres,
           'p_thres|p=i' => \$n_thres,
           'outfile|o=s' => \$outfile,
           'mode|m=s' => \$mode,
           'help|h' => \my $help,
           'verbose|v' => \my $verbose);

help() if ( defined $help );
verbose() if(defined $verbose or !defined($file) or !defined($start_ind) or !defined($end_ind));
print STDERR "Options USED:
file(f)=$file\toutfile(o)=$outfile\tmode(m)=$mode
chr(c)=$chr\t\tstart(s)=$start\t\tend(e)=$end
start_ind(a)=$start_ind\t\tend_ind(b)=$end_ind\t\tc_min(n)=$control_min\t\tc_max(x)=$control_max
c_q1(1)=$control_q1\tc_q9(9)=$control_q9\tc_mean(m)=$control_mean\tc_mean_thres(u)=$control_mean_thres
beta_thres(t)=$b_thres\tp_thres(p)=$n_thres\t\twindow(w)=$window\n";
if(($mode ne "BOTH") && ($mode ne "MAX_ONLY") && ($mode ne "MIN_ONLY")){ print STDERR "Warning!! Invalid option $mode. Overriding with BOTH\n"}

#open(IN, $file) or die "***Error*** Can not open the file $file";
my $i = 0; my $j = 0;
#my $header = <IN>;
my @data = ();
my $is_stdin = 0;
my $IN;
if($file eq "stdin" || $file eq "STDIN"){
    $IN = *STDIN;
    $is_stdin++;
}else{
    open $IN, "<", $file or die $!
}

my $header = <$IN>; chomp $header; my @h =  split(/\t/,$header);
while(<$IN>){
    chomp $_;
    my @s = split(/\t/,$_);
    for($j = 0; $j<@s; $j++){
        $data[$i][$j] = $s[$j];
    }
    $i++;
}
close $IN unless $is_stdin;

my $nrow = $i;
my $ncol = $j;
undef $i;
undef $j;

my @sign = (0) x $nrow;
my @sign_ind = ('NA') x $nrow;
my @mean = ('NA') x $nrow;
my @direction_all = ('NA') x $nrow;
my @num = (0) x $nrow;
my $corr = 0;#added for correction at the end of each chromosome
my @xx = ('NA') x $nrow;

for(my $k=$start_ind; $k<=$end_ind; $k++){
    my @direction = ('NA') x $nrow;
    for(my $i = 0; $i<$nrow; $i++){
        my $count = 0; my $count1=0; my $count2=0; my $total_p = 0; my $sum = 0;
        my $countq1=0; my $countq2=0;
        for(my $j = $i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$j+1][$chr]; $j++){
#       print "$j P = $data[$j][$p] \n";
            if(looks_like_number($data[$j][$k])){
                if(($data[$j][$k] >= $data[$j][$control_q9]) && ($data[$j][$k] >= $data[$j][$control_mean] + $control_mean_thres)){ $countq1++; }
                if(($data[$j][$k] <= $data[$j][$control_q1]) && ($data[$j][$k] <= $data[$j][$control_mean] - $control_mean_thres)){ $countq2++; }
                if($data[$j][$k] >= $data[$j][$control_max] + $b_thres){ $count1++; }
                if($data[$j][$k] <= $data[$j][$control_min] - $b_thres){ $count2++; }
                $total_p++;
                $corr = $j;
            }
        }

    #Correction for the change in chromosome number
        if($total_p==0){  $corr = $i+1;}
        if($total_p==0){
#       print $j,"\t",$corr,"\t",$data[$corr-1][$chr],"\t", $data[$corr][$chr],"\n";
            if(looks_like_number($data[$j][$k])){
                if(($data[$corr-1][$k] >= $data[$corr-1][$control_q9]) && ($data[$j][$k] >= $data[$j][$control_mean] + $control_mean_thres)){ $countq1++; }
                if(($data[$corr-1][$k] <= $data[$corr-1][$control_q1]) && ($data[$j][$k] <= $data[$j][$control_mean] - $control_mean_thres)){ $countq2++; }
                if($data[$corr-1][$k] >= $data[$corr-1][$control_max] + $b_thres){ $count1++; }
                if($data[$corr-1][$k] <= $data[$corr-1][$control_min] - $b_thres){ $count2++; }
                $total_p++;
            }
        }

        my $flag='';
        if($mode eq "BOTH"){ $flag = ((($countq1>=$n_thres) && ($count1>=1))||(($countq2>=$n_thres) && ($count2>=1))); }
        elsif($mode eq "MAX_ONLY") {$flag = (($countq1>=$n_thres) && ($count1>=1)); }
        elsif($mode eq "MIN_ONLY") {$flag = (($countq2>=$n_thres) && ($count2>=1)); }
        else{ $flag = ((($countq1>=$n_thres) && ($count1>=1))||(($countq2>=$n_thres) && ($count2>=1))); }

        if($flag){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$j+1][$chr]; $j++){
                $sign[$j] = 1;
                $h[$k] =~ s/.AVG_Beta$//;
                my $temp = ",",$h[$k].",";
                my @tt = split(",", $sign_ind[$j]);
                my $match_flag = 0;
                if($sign_ind[$j] eq "NA") {$sign_ind[$j] = $h[$k];}
                elsif(!($h[$k]~~@tt)){$sign_ind[$j] .= ",".$h[$k];}
#               elsif(($sign_ind[$j] !~ /^$h[$k]/) && ($sign_ind[$j] !~ /$temp/) && ($sign_ind[$j] !~ /$h[$k]$/)) {$sign_ind[$j] .= ",".$h[$k];}
            }
        }
        if((($countq1>=$n_thres) && ($count1>=1))&&(($countq2>=$n_thres) && ($count2>=1))){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$j+1][$chr]; $j++){ $direction[$j] ="BOTH";}
        }elsif((($countq1>=$n_thres) && ($count1>=1)) && ($mode eq "BOTH" || $mode eq "MAX_ONLY")){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$j+1][$chr]; $j++){ $direction[$j] ="MAX";}
        }elsif((($countq2>=$n_thres) && ($count2>=1)) && ($mode eq "BOTH" || $mode eq "MIN_ONLY")){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$j+1][$chr]; $j++){ $direction[$j] ="MIN";}
        }
    }
    for(my $i = 0; $i<@direction; $i++){
        if($direction_all[$i] eq "NA"){ $direction_all[$i] = $direction[$i]}
        else{
            if($direction[$i] ne "NA"){ $direction_all[$i] .= ",".$direction[$i]}
        }
    }

}

my $OUT;
my $is_stdout = 0;
if($outfile eq "stdout" || $outfile eq "STDOUT"){
    $is_stdout++;
}else{
    open(STDOUT, ">".$outfile);
}

#perl window2.pl aa 1 2 3 1000 0 18 1 0.1 3 0.5
print STDOUT $header,"\tSign_individuals_t",$b_thres,"_n",$n_thres,"_w",$window/1000,"k_",$mode,"\tSign_direction_t",$b_thres,"_n",$n_thres,"_w",$window/1000,"k_",$mode,"\tSign_window_t",$b_thres,"_n",$n_thres,"_w",$window/1000,"k_",$mode,"\n";

for(my $i = 0; $i <$nrow; $i++){
    for(my $j = 0; $j<$ncol; $j++){
        print STDOUT $data[$i][$j],"\t";
    }
    print STDOUT $sign_ind[$i],"\t",$direction_all[$i],"\t",$sign[$i],"\n";
}

close STDOUT unless $is_stdout;

sub help{
    print STDERR "Usage: perl .pl -f filename -a 30 -b 122 -o output
perl window_analysis.pl -f filename -c chr_col -s start_col -e end_col -w window -a ind_start -b ind_end  -n control_min -x control_max -t beta_thres -p probe_thres -o outfile -m mode\n";
    exit();
}

sub verbose{
print STDERR "Usage: perl .pl -f filename -a 30 -b 122 -o output
perl window_analysis.pl -f filename -c chr_col -s start_col -e end_col -w window -a ind_start -b ind_end  -n control_min -x control_max -t beta_thres -p probe_thres -ooutfile -m mode

-f, --filename=FILENAME\t\t\tInput filename sorted  by chromosome and probe start, can be stdin or STDIN while piping (REQUIRED)
-a, --start_ind=IND_START_COL\t\tColumn from where sample beta value starts (REQUIRED)
-b, --end_ind=IND_END_COL\t\tColumn from where sample beta value ends (REQUIRED)
-o, --outfile=OUTFILE\t\t\tOutput file name with three extra columns for Significant indiduals, directionality and window,
\t\t\t\t\tcan be STDOUT or stdout for printing commandline, can be used for piping and can be redirected to a file using '>'
\t\t\t\t\t(default = STDOUT)
-c, --chr=CHR_COL\t\t\tColumn containing probe chromosome information (default = 1)
-s, --start=START_COL\t\t\tColumn containing probe start information (default = 2)
-e, --end=END_COL\t\t\tColumn containing probe end information (default = 3)
-w, --window=WINDOW_SIZE\t\tSize of the sliding window (default = 1000)
-n, --c_min=CONTROL_MIN_COL\t\tColumn which contains control min beta values (default = 24)
-x, --c_max=CONTROL_MIN_COL\t\tColumn which contains control max beta values (default = 23)
-t, --beta_thres=BETA_THRESHOLD\t\tBeta value threshold +/- of control min/max (default = 0.1)
-p, --p_thres=NO_OF_PROBES_THRESHOLD\tNumber of significant probes in window threshold (default = 2)
-m, --mode=MODE\t\t\t\tCan be BOTH, MAX_ONLY or MIN_ONLY. runs script in both or one direction (default=BOTH)
-h, --help\t\t\t\tUsage summary
-v, --verbose\t\t\t\tDetailed Usage Information\n";
exit();
}
