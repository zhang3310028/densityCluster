#!/usr/bin/perl -w
use Math::Round qw/round/;
use List::Util;
use strict;
use Getopt::Long;
use Statistics::R;
use Time::HiRes qw(time);
### how to decide rho and delta?
my ($data_file,$outp,$t_value,$k_value,$verbose,$delta_rho,$gaussian_kenel,$has_header,$cols,$rho_min,$delta_min,$cols_file) =  
("","cluster.out",0.02,-1,"","","","","",0,0,"") ;
GetOptions ( 
	"t_value=f" => \$t_value,    # float
	"k=i" => \$k_value,
	"o=s" => \$outp,
	"g" => \$gaussian_kenel,
	"H" => \$has_header,
	"file=s"   => \$data_file,      # string
	"o_rho_file=s" => \$delta_rho,
	"cols_file=s" => \$cols_file,
	"cols=s"=> \$cols,
	"rho_min=f" => \$rho_min,
	"delta_min=i" => \$delta_min,
	"verbose"  => \$verbose)   # flag
or die("Error in command line arguments\n");


my $cols_data_ref;
my $data_num ;
my $cutoff_distance;
my $cols_ref;
my $cur_time = time();
my @col_lines; # 代表每行数据的hotspot位置
sub init{
	if($verbose){
		$cur_time = time() ;
		warn "\ninit ..."."\n";
	}
	if($cols_file){
		open(COLSFILE,"<$cols_file");
		@col_lines = <COLSFILE>;
		chomp(@col_lines);
		close(COLSFILE);
	}else{
		$cols !~ /(\d+)?(,\d+)+?/ and die ("Error: parameter -cols");
	}
	
	open(F,"<$data_file");
	my $header ;
	if($has_header){
		$header = <F>;
		chomp $header;
	}
	my @data = <F>;
	close(F);

	$data_num = scalar(@data);
	warn "data_number:$data_num"."\n" if($verbose);

	chomp(@data);
	if($cols_file){
		$cols_data_ref = &get_data_cols_exactly(\@data,\@col_lines);
	}else{
		my @cols;
		if($cols eq ""){
			my $tmp_line = $data[0];
		}else{
			@cols = split(/,/,$cols);
		}
		$cols_data_ref = &get_data_cols(\@data,\@cols);
		$cols_ref = \@cols;
	}
	if($verbose){
		my $timestamp = time();
		my $time_used = $timestamp - $cur_time;
		$cur_time = $timestamp;
		warn "init time used:$time_used\n";
	}
}
sub get_data_cols{
	my ($data_ref,$cols_ref)=@_;
	my $cols_data_ref=[] ;
	for(my $i = 0 ; $i < @{$data_ref} ;$i ++){
		my $line = $data_ref->[$i];
		$line  =~ s/\s*$//g;
		my @tmp = split(/\s+/,$line);
		my @tmp_col = @tmp[@{$cols_ref}];
		$cols_data_ref->[$i] = \@tmp_col;
	}
	return $cols_data_ref;
}
sub get_data_cols_exactly(){
	my ($data_aref,$col_lines_aref)=@_;
	my $cols_data_ref=[] ;
	for(my $i = 0 ; $i < @{$data_aref} ;$i ++){
		my $line = $data_aref->[$i];
		my $cols = $col_lines_aref->[$i];
		my @cols=split(/\t/,$cols);
		my @cols_index;
		foreach $i(0 .. $#cols){
			$cols_index[$i] = $cols[$i]-1;
		}
		$line  =~ s/\s*$//g;
		my @tmp = split(/\s+/,$line);
		my @tmp_col = @tmp[@cols_index];
		$cols_data_ref->[$i] = \@tmp_col;
	}
	return $cols_data_ref;
}

my %delta_ends=();
&main();
sub main{
	&init();
	if($verbose){
		warn  "step1 findClustrMedoids_index() ...\n";
	}
#	my ($cluster_medoids_index_ref,$rho_ref,$pdistances_aref,$index_by_desc_rho_ref) = &findClustrMedoids_index();
	my ($cluster_medoids_index_ref,$rho_ref,$index_by_desc_rho_ref) = &findClustrMedoids_index();
	if($verbose){
		warn  "classifycluster() ...\n";
	}
	my $cluster_ref = &Classifycluster($cluster_medoids_index_ref,$index_by_desc_rho_ref);
	if($verbose){
		warn  "decide_clusterCore() ...\n";
	}
#	my $halo_ref=&decide_clusterCore($cluster_ref,$rho_ref,$pdistances_aref);
	my $halo_ref=&decide_clusterCore2($cluster_ref,$rho_ref,$cols_data_ref);

	open(F,">$outp.cluster");
	print F join("\n",@{$cluster_ref});
	close(F);

	open(F,">$outp.halo.txt");
	print F join("\n",@{$halo_ref});
	close(F);

	open(F,">$outp.cluster_medoids.txt");
	print F join("\n",@{$cluster_medoids_index_ref});
	close(F);
}

sub decide_clusterCore(){
	my ($cluster_ref,$rho_ref,$pdistances_aref) = @_;
	my @halo;
	my @border_rho;
	for my $index (0..(@{$cluster_ref}-1)){
		$border_rho[$index] = 0;
	}
	for my $key1(keys @{$pdistances_aref}){
		for my $key2 (keys @{$pdistances_aref->[$key1]}){
			next if (!defined $pdistances_aref->[$key1]->[$key2]);
			if($pdistances_aref->[$key1]->[$key2] < $cutoff_distance && $cluster_ref->[$key1] != $cluster_ref->[$key2]){
				my $rho_aver = ($rho_ref->{$key1} + $rho_ref->{$key2})/2;
				if($rho_aver > $border_rho[$cluster_ref->[$key1]]){
					$border_rho[$cluster_ref->[$key1]] = $rho_aver;
				}
				if($rho_aver > $border_rho[$cluster_ref->[$key2]]){
					$border_rho[$cluster_ref->[$key2]] = $rho_aver;
				}

			}
		}
	}
	for(my $i = 0 ; $i < $data_num ; $i ++){
		if($rho_ref->{$i} < $border_rho[$cluster_ref->[$i]]){
			$halo[$i] = 0;
		}else{
			$halo[$i] = 1;
		}
	}
	return \@halo;
}

sub decide_clusterCore2(){
	my ($cluster_ref,$rho_ref,$cols_data_aref) = @_;
	my @halo;
	my @border_rho;
	for my $index (0..(@{$cluster_ref}-1)){
		$border_rho[$index] = 0;
	}
	my $data_num = @{$cols_data_aref};
	for(my $i = 0 ; $i < $data_num - 1 ; $i ++){
		for(my $j = $i+1 ; $j < $data_num ; $j ++){
			if(&calculate_ppdistance2($cols_data_aref,$i,$j) < $cutoff_distance && $cluster_ref->[$i] != $cluster_ref->[$j]){
				my $rho_aver = ($rho_ref->{$i} + $rho_ref->{$j})/2;
				if($rho_aver > $border_rho[$cluster_ref->[$i]]){
					$border_rho[$cluster_ref->[$i]] = $rho_aver;
				}
				if($rho_aver > $border_rho[$cluster_ref->[$j]]){
					$border_rho[$cluster_ref->[$j]] = $rho_aver;
				}
			}
		}
	}

	for(my $i = 0 ; $i < $data_num ; $i ++){
		if($rho_ref->{$i} < $border_rho[$cluster_ref->[$i]]){
			$halo[$i] = 0;
		}else{
			$halo[$i] = 1;
		}
	}
	return \@halo;
}

sub findClustrMedoids_index(){
	if($verbose){
		warn "step1.1 calculate_pdistance() ...\n";
		$cur_time = time();
	}
#	my $pdistances_aref = &calculate_pdistance($cols_data_ref);
	if($verbose){
		my $timestamp = time();
		my $time_used = $timestamp - $cur_time;
		$cur_time = $timestamp;
		warn "time used:$time_used\n";
		warn "step1.2 calculate_cutoff() ...\n";
	}
#	$cutoff_distance = &calculate_cutoff($pdistances_aref,$t_value);
#	$cutoff_distance = &calculate_cutoff2($cols_data_ref,$t_value);
#	$cutoff_distance = &calculate_cutoff2_1($cols_data_ref,$t_value);
	$cutoff_distance = &calculate_cutoff2_2($cols_data_ref,$t_value);
	if($verbose){
		my $timestamp = time();
		my $time_used = $timestamp - $cur_time;
		$cur_time = $timestamp;
		warn "time used:$time_used\n";
		warn "cutoff_distance = $cutoff_distance...\n";
	}
	open(F,">$outp.cutoff_distance.txt");
	print F $cutoff_distance;
	close(F);
	if($verbose){
		my $timestamp = time();
		my $time_used = $timestamp - $cur_time;
		$cur_time = $timestamp;
		warn "time used:$time_used\n";
		warn "step1.3 calculate_rho() ...\n";
	}
#	my $rho_ref = &calculate_rho($pdistances_aref,$cutoff_distance);
#	my $rho_ref = &calculate_rho2($cols_data_ref,$cutoff_distance);
#	my $rho_ref = &calculate_rho2_1($cols_data_ref,$cutoff_distance);
#	my $rho_ref = &calculate_rho3_1($cols_data_ref,$cutoff_distance);
	my $rho_ref = &calculate_rho3_2($cols_data_ref,$cutoff_distance);
	if($verbose){
		my $timestamp = time();
		my $time_used = $timestamp - $cur_time;
		$cur_time = $timestamp;
		warn "time used:$time_used\n";
		warn "step1.4 get_index_by_desc_rho() ...\n";
	}
	my $index_by_desc_rho_ref = &get_index_by_desc_rho($rho_ref);
	my $deltas_ref = &calculate_delta($index_by_desc_rho_ref,$rho_ref);
	my $index_of_deltaRho_desc_ref=&calc_index_of_deltaRho_desc($rho_ref,$deltas_ref);
	my @cluster_medoids_index;
	if($k_value != -1){
		my $cur_cluster = 0;
		$cluster_medoids_index[0] = $index_by_desc_rho_ref->[0];
		for(my $cur_index = 1 ; $cur_index < @{$index_by_desc_rho_ref} ; $cur_index ++){
			if($rho_ref->{$index_by_desc_rho_ref->[$cur_index]} >= $rho_min 
				&& $deltas_ref->{$index_by_desc_rho_ref->[$cur_index]} >= $delta_min){
				$cluster_medoids_index[++$cur_cluster] = $index_by_desc_rho_ref->[$cur_index];
				$cur_cluster+1 >=$k_value  and last;
			}
		}
			
	}else{
		my $cur_cluster = 1;
		$cluster_medoids_index[0] = $index_by_desc_rho_ref->[0];
		my $cur_delta = 0;
		
		for(my $i = 1 ; $i < @{$index_by_desc_rho_ref} ; $i ++){
			if($cur_delta!=0 && $deltas_ref->{$index_by_desc_rho_ref->[$i]} / ($cur_delta/($i-$cur_cluster)) > 30){
				$cluster_medoids_index[$cur_cluster++] = $index_by_desc_rho_ref->[$i];
			}else{
				$cur_delta += $deltas_ref->{$index_by_desc_rho_ref->[$i]};
			}
		}
	
	}
#	return (\@cluster_medoids_index,$rho_ref,$pdistances_aref,$index_by_desc_rho_ref);
	return (\@cluster_medoids_index,$rho_ref,$index_by_desc_rho_ref);
}

sub get_index_by_desc_rho(){
	my $rho_ref = shift;
	my @indexs_by_desc_rho;
	@indexs_by_desc_rho = sort{$rho_ref->{$b} <=> $rho_ref->{$a}} keys %{$rho_ref};
	return \@indexs_by_desc_rho;
}

sub calc_index_of_deltaRho_desc{
	my ($rho_ref,$deltas_ref) = @_;
	my @index_result;
	@index_result = sort{
		$rho_ref->{$b} * ($deltas_ref->{$b}) <=> $rho_ref->{$a} * ($deltas_ref->{$a}) 
#		($rho_ref->{$b}) <=> ($rho_ref->{$a}) 
#		($deltas_ref->{$b}) <=> ($deltas_ref->{$a}) 
	} keys %{$deltas_ref};
	return \@index_result;
}

sub calculate_delta{
	my ($index_by_rho_desc_ref,$rho_ref) = @_;
	my %deltas;
	my $max_dis = -1;
	for my $i (0 .. ($data_num-1)){
		$deltas{$i} = 0;
	}
	
	for(my $i = 1 ; $i < @{$index_by_rho_desc_ref} ; $i ++){
		my $cur_delta = -1;
		my $cur_delta_end = -1;
		for(my $j = 0 ; $j < $i ; $j ++){
			my $distance = &calculate_ppdistance($index_by_rho_desc_ref->[$i],$index_by_rho_desc_ref->[$j]);
			if($cur_delta == -1 || $distance < $cur_delta){
				$cur_delta = $distance;
				$cur_delta_end = $index_by_rho_desc_ref->[$j];
			}
		}
		$delta_ends{$index_by_rho_desc_ref->[$i]} = $cur_delta_end;
		$deltas{$index_by_rho_desc_ref->[$i]} = $cur_delta;
		if($cur_delta > $max_dis){
			$max_dis = $cur_delta;
		}

	}
	$delta_ends{$index_by_rho_desc_ref->[0]} = $index_by_rho_desc_ref->[0];
	$deltas{$index_by_rho_desc_ref->[0]} = $max_dis;

	if($verbose){
		my $count = scalar(keys %deltas);
		warn "delta_count:".$count."\n"
	}
	open(F,">$outp\_rho_delta.txt");
	for my $key (sort {$a <=> $b} keys %{$rho_ref}){
		print F "$key\t".$rho_ref->{$key}."\t$deltas{$key}"."\n";
	}
	close(F);
	return \%deltas;
}

sub calculate_rho{
	my ($pdistances_aref,$cutoff_distance) = @_ ;
	my %rho_values;
	for (my $i=0 ; $i <$data_num ; $i++){
		$rho_values{$i} = 0;
	}
	foreach my $point_s (keys @{$pdistances_aref}){
		foreach my $point_e (keys @{$pdistances_aref->[$point_s]}){
			my $distance = $pdistances_aref->[$point_s]->[$point_e];
			next if (!defined $distance);
			if($gaussian_kenel){
				$rho_values{$point_s} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
				$rho_values{$point_e} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
			}else{
				if($distance < $cutoff_distance){
					$rho_values{$point_s} ++;
					$rho_values{$point_e} ++;
				}
			}

			
		}
	}
	if($verbose){
		my $rho_count = scalar(keys %rho_values);
		print "rho_count:" . $rho_count."\n";
	}
	return \%rho_values;
}

sub calculate_rho2{ # 避免大距离矩阵的计算计算距离
	my ($col_data_aref,$cutoff_distance) = @_ ;
	my %rho_values;
	my $data_num = @{$col_data_aref};
	for (my $i=0 ; $i <$data_num ; $i++){
		$rho_values{$i} = 0;
	}
	for(my $i = 0 ; $i < $data_num - 1 ; $i ++){
		for(my $j = $i+1 ; $j < $data_num ; $j ++){
			my $distance = &calculate_ppdistance2($col_data_aref,$i,$j);
			if($gaussian_kenel){
				$rho_values{$i} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
				$rho_values{$j} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
			}else{
				if($distance <= $cutoff_distance){
					if($distance < $cutoff_distance){
						$rho_values{$i} ++;
						$rho_values{$j} ++;
					}
				}
			}
		}
	}
	if($verbose){
		my $rho_count = scalar(keys %rho_values);
		print "rho_count:" . $rho_count."\n";
	}
	return \%rho_values;
}

sub calculate_rho2_1{
	my ($col_data_aref,$cutoff_distance) = @_ ;
	my %rho_values;
	my $data_num = @{$col_data_aref};
	for (my $i=0 ; $i <$data_num ; $i++){
		$rho_values{$i} = 0;
	}
	for(my $i = 0 ; $i < $data_num - 1 ; $i ++){
		for(my $j = $i+1 ; $j < $data_num ; $j ++){
			my $distance = &calculate_ppdistance2($col_data_aref,$i,$j);
			if($distance <= $cutoff_distance){
				if($gaussian_kenel){
					$rho_values{$i} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
					$rho_values{$j} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
				}else{
						$rho_values{$i} ++;
						$rho_values{$j} ++;
				}
			}
		}
	}
	if($verbose){
		my $rho_count = scalar(keys %rho_values);
		print "rho_count:" . $rho_count."\n";
	}
	return \%rho_values;
}
sub calculate_rho3_1{
	my ($col_data_aref,$cutoff_distance) = @_ ;
	my %rho_values;
	my $data_num = @{$col_data_aref};
	my $R = Statistics::R->new();
	for (my $i=0 ; $i <$data_num ; $i++){
		$rho_values{$i} = 0;
	}
	my ($sample_cur,$sample_count) = (0,$data_num * 100);
	$R->startR();
	my $data_end_index = $data_num - 1;
	while($sample_cur < $sample_count){
		$R->send(qq`sample_index<-sample(0:$data_end_index,2)`);
		my $index_aref=$R->get(qq`sample_index`);
		my $distance = &calculate_ppdistance2($col_data_aref,$index_aref->[0],$index_aref->[1]);
		if($distance <= $cutoff_distance){
				if($gaussian_kenel){
					$rho_values{$index_aref->[0]} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
					$rho_values{$index_aref->[1]} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
				}else{
						$rho_values{$index_aref->[0]} ++;
						$rho_values{$index_aref->[1]} ++;
				}
		}
		$sample_cur++;
	}
	
	if($verbose){
		my $rho_count = scalar(keys %rho_values);
		print "rho_count:" . $rho_count."\n";
	}
	$R->stopR();
	return \%rho_values;
}

sub calculate_rho3_2{
	my ($col_data_aref,$cutoff_distance) = @_ ;
	my %rho_values;
	my $data_num = @{$col_data_aref};
	for (my $i=0 ; $i <$data_num ; $i++){
		$rho_values{$i} = 0;
	}
	my ($sample_cur,$sample_count) = (0,$data_num * 100);
	while($sample_cur < $sample_count){
		my @index = &sample_no_rep($data_num,2);
		my $distance = &calculate_ppdistance2($col_data_aref,$index[0],$index[1]);
		if($distance <= $cutoff_distance){
				if($gaussian_kenel){
					$rho_values{$index[0]} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
					$rho_values{$index[1]} += exp(-($distance/$cutoff_distance)*($distance/$cutoff_distance));
				}else{
						$rho_values{$index[0]} ++;
						$rho_values{$index[1]} ++;
				}
		}
		$sample_cur++;
	}
	
	if($verbose){
		my $rho_count = scalar(keys %rho_values);
		print "rho_count:" . $rho_count."\n";
	}
	return \%rho_values;
}

sub Classifycluster(){
	my ($cluster_medoids_index_ref,$index_by_desc_rho_ref) = @_;
	my @cluster;
	for(my $i = 0 ; $i < @{$cluster_medoids_index_ref} ;$i ++){
		$cluster[$cluster_medoids_index_ref->[$i]]=$i;
		$delta_ends{$cluster_medoids_index_ref->[$i]}=$cluster_medoids_index_ref->[$i];
	}

	for(my $i = 0 ; $i < @{$index_by_desc_rho_ref} ; $i++){
		my $end_index = $delta_ends{$index_by_desc_rho_ref->[$i]};
		my $cur_cluster =$cluster[$end_index];
		$cluster[$index_by_desc_rho_ref->[$i]] = $cur_cluster;
	}
	return \@cluster;
}
sub calculate_pdistance{
	my $pdistance_aref=[];
	my $data_num = @{$cols_data_ref};
	for(my $i = 0 ; $i < $data_num-1 ; $i ++){
		$pdistance_aref->[$i] = [];
		for(my $j = $i+1 ; $j < $data_num; $j ++){
			my $dist = &calculate_ppdistance($i,$j);
			$pdistance_aref->[$i]->[$j] = $dist;
		}
	}
	return $pdistance_aref;
}

sub calculate_ppdistance{
	my ($data_index_1,$data_index_2) = @_;
	my $line1_ref = $cols_data_ref->[$data_index_1];
	my $line2_ref = $cols_data_ref->[$data_index_2];
	my $distance = 0;
	for(my $i = 0 ; $i < @{$cols_data_ref->[0]} ; $i ++){
		if(not defined $line2_ref->[$i]){
			$line2_ref->[$i] = 0;
		}
		if(not defined $line1_ref->[$i]){
			$line1_ref->[$i] = 0;
		}
		$distance += ($line1_ref->[$i]-$line2_ref->[$i])**2;
	}
	$distance = sqrt($distance);
	return $distance;
}
sub calculate_ppdistance2{
	my ($cols_data_aref,$data_index_1,$data_index_2) = @_;
	my $line1_ref = $cols_data_aref->[$data_index_1];
	my $line2_ref = $cols_data_aref->[$data_index_2];
	my $distance = 0;
	for(my $i = 0 ; $i < @{$cols_data_aref->[0]} ; $i ++){
		if(not defined $line2_ref->[$i]){
			$line2_ref->[$i] = 0;
		}
		if(not defined $line1_ref->[$i]){
			$line1_ref->[$i] = 0;
		}
		$distance += ($line1_ref->[$i]-$line2_ref->[$i])**2;
	}
	$distance = sqrt($distance);
	return $distance;
}
sub calculate_cutoff{
	my ($pdistances_aref , $t_value) = @_;
	my @sorted_distance = &sort_distance($pdistances_aref);
	my $distance_num = scalar(@sorted_distance);
	if($distance_num != $data_num*($data_num-1)/2){
		print STDERR "data_num:".$data_num.",distance_num:" .$distance_num;
		print STDERR "distance not match data number!" and die;
	}
	my $k = &round($distance_num * $t_value);
	my $cutoff_distance = $sorted_distance[$k-1];
	if($cutoff_distance == 0){
		$cutoff_distance = 0.0000000001;
	}
	return $cutoff_distance;
}
sub calculate_cutoff2{
	my ($col_data_aref , $t_value) = @_;
	my $k = &round($data_num*($data_num-1)/2 * $t_value);
#	my @sorted_distance = &sort_distance1($col_data_aref);
	my @sorted_distance = &sort_distance2($col_data_aref,$k);
	my $distance_num = scalar(@sorted_distance);
	my $cutoff_distance = $sorted_distance[$k-1];
	if($cutoff_distance == 0){
		$cutoff_distance = 0.0000000001;
	}
	return $cutoff_distance;
}

sub sort_distance{
	my ($pdistances_aref) = @_ ;
	my @sorted_distance ;
	for my $k1 (keys @{$pdistances_aref}){
		for my $k2(keys @{$pdistances_aref->[$k1]}){
			! defined $pdistances_aref->[$k1]->[$k2] and next;
			$sorted_distance[@sorted_distance]= $pdistances_aref->[$k1]->[$k2];
		}
	}
	@sorted_distance = sort{$a <=> $b} @sorted_distance;
	return @sorted_distance;
}
sub sort_distance1{
	my ($col_data_aref) = @_ ;
	my @sorted_distance ;
	my $data_num = @{$col_data_aref};
	for(my $i = 0 ; $i < $data_num - 1 ; $i ++){
		for(my $j = $i+1 ; $j < $data_num ; $j ++){
			my $dist= &calculate_ppdistance2($col_data_aref,$i,$j);
			$sorted_distance[@sorted_distance]=  $dist;
		}
	}
	@sorted_distance = sort{$a <=> $b} @sorted_distance;
	return @sorted_distance;
}
sub sort_distance2{
	my ($col_data_aref,$kvalue) = @_ ;
	my @sorted_distance;
	my $data_num = @{$col_data_aref};
	my ($cur_index,$cur_max_dist) = (-1,-1);
	my $last_arr_len = 0;
	for(my $i = 0 ; $i < $data_num - 1 ; $i ++){
		for(my $j = $i+1 ; $j < $data_num ; $j ++){
			my $dist= &calculate_ppdistance2($col_data_aref,$i,$j);
			my $dis_num = @sorted_distance;
			if($dis_num < $kvalue){
				$last_arr_len = $dis_num;
				$sorted_distance[$dis_num] = $dist;
				if($dis_num +1 == $kvalue){
					@sorted_distance = sort{$a <=> $b} @sorted_distance;
					$cur_index = @sorted_distance-1;
					$cur_max_dist = $sorted_distance[$cur_index];
				}
			}elsif($dist < $cur_max_dist){
				my $insert_index = &get_insert_index($dist,@sorted_distance);
				splice(@sorted_distance,$insert_index,0,$dist);
				$#sorted_distance =  $kvalue;
				$cur_max_dist = $sorted_distance[$cur_index];
			}

		}
	}
	@sorted_distance = sort{$a <=> $b} @sorted_distance;
	return @sorted_distance;
}
sub calculate_cutoff2_1{ # 取样解决计算cutoff
	my ($col_data_aref,$tvalue) = @_ ;
	my $R = Statistics::R->new();
	$R->startR();
	my @sorted_distance;
	my $data_num = @{$col_data_aref};
	my $dist_num = $data_num*($data_num-1)/2;
	my ($sample_cur, $sample_count) = (0,1000);
	if($sample_count > $dist_num){ #TODO 这一步可以不用取样
		$sample_count = $dist_num;
	}
	my $kvalue = $dist_num * $t_value;
	my ($mdist_index,$cur_mdist) = (-1,-1);
	my $last_arr_len = 0;
	my $data_end_index= $data_num - 1;
	my $sample_kvalue = &round($kvalue * $sample_count / $dist_num); #TODO sample_kvalue == 0 
#	die " $dist_num,$kvalue,$sample_count,$sample_kvalue";
	my @sample_sorted_distance=();
	while($sample_cur < $sample_count){
		$R->send(qq`sample_index<-sample(0:$data_end_index,2)`);
		my $index_aref=$R->get(qq`sample_index`);
		#	die join(",",@{$index_aref});
		my $dist = &calculate_ppdistance2($col_data_aref ,$index_aref->[0], $index_aref->[1]);
		$sample_sorted_distance[@sample_sorted_distance] = $dist;
		$sample_cur ++ ;
	}
	@sample_sorted_distance = sort{$a <=> $b} @sample_sorted_distance;
	my $cutoff_distance = $sample_sorted_distance[$sample_kvalue-1];
	$R->stopR();
	return $cutoff_distance;
}
sub calculate_cutoff2_2{ # 避免R抽样
	my ($col_data_aref,$tvalue) = @_ ;
	my @sorted_distance;
	my $data_num = @{$col_data_aref};
	my $dist_num = $data_num*($data_num-1)/2;
	my ($sample_cur, $sample_count) = (0,1000);
	if($sample_count > $dist_num){ #TODO 这一步可以不用取样
		$sample_count = $dist_num;
	}
	my $kvalue = $dist_num * $t_value;
	my ($mdist_index,$cur_mdist) = (-1,-1);
	my $last_arr_len = 0;
	my $data_end_index= $data_num - 1;
	my $sample_kvalue = &round($kvalue * $sample_count / $dist_num); #TODO sample_kvalue == 0 
#	die " $dist_num,$kvalue,$sample_count,$sample_kvalue";
	my @sample_sorted_distance=();
	while($sample_cur < $sample_count){
		my @index = &sample_no_rep($data_num,2);
		my $dist = &calculate_ppdistance2($col_data_aref ,$index[0], $index[1]);
		$sample_sorted_distance[@sample_sorted_distance] = $dist;
		$sample_cur ++ ;
	}
	@sample_sorted_distance = sort{$a <=> $b} @sample_sorted_distance;
	my $cutoff_distance = $sample_sorted_distance[$sample_kvalue-1];
	return $cutoff_distance;
}

sub get_max_index{

	my @arr = @_;
	my $index = 0;
	my $max_value = $arr[0];
        foreach(1 .. $#arr)
        {
            if($arr[$_] > $max_value)
            {
                $index = $_;
		$max_value = $arr[$_];
            }
        }
	return $index;
}
sub get_insert_index{
	my ($dist,@sorted_distance_asc) = @_;
	my $num_dis = @sorted_distance_asc;
	my $index = $num_dis-1;
	while($index >= 0 && $sorted_distance_asc[$index] > $dist){
		$index -- ;
	}
	return $index+1;
}
sub sample_no_rep {
	my ($n, $m) = @_;
	return if ($m > $n);
	my @b;
	if ($m > sqrt($n)) { 
		my $i = $n;
		for (my $j = $m; $j > 0; --$j) {
			my ($p, $x) = (1.0, rand());
			$p -= $p * $j / ($i--) while ($x < $p);
			push(@b, $n - $i - 1);
		}
		return @b;
	} else { # small $m
		my %s;
		while (@b < $m) {
			$_ = int(rand($n));
			unless ($s{$_}) {
				push(@b, $_);
				$s{$_} = 1;
			}
		}
		return @b;
	}
}
