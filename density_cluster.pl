#!/usr/bin/perl -w
use Math::Round qw/round/;
use strict;
use Getopt::Long;
### how to decide rho and delta?
my ($data_file,$out_file,$t_value,$k_value,$verbose,$delta_rho,$gaussian_kenel,$has_header)=("","cluster.out",0.02,-1,"","","","") ;
GetOptions ( "t_value=f" => \$t_value,    # float
	"k=i" => \$k_value,
	"o=s" => \$out_file,
	"g" => \$gaussian_kenel,
	"H" => \$has_header,
	"file=s"   => \$data_file,      # string
	"o_rho_file=s" => \$delta_rho,
	"verbose"  => \$verbose)   # flag
or die("Error in command line arguments\n");

my $cols_data_ref;
my $data_num ;
my $cutoff_distance;
my $cols_ref;
sub init(){
	open(F,"<$data_file");
	my $header ;
	if($has_header){
		$header = <F>;
		chomp $header;
	}
	my @data = <F>;
	my @cols = (121,124);
	$cols_ref = \@cols;
	chomp(@data);
	close(F);
	$cols_data_ref = &get_data_cols(\@data,\@cols);
	$data_num = scalar(@data);
}
sub get_data_cols(){
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

my %delta_ends=();
&main();
sub main(){
	&init();
	if($verbose){
		print  "step1 findClustrMedoids_index() ...\n";
	}
	my ($cluster_medoids_index_ref,$rho_ref,$pdistances_ref,$index_by_desc_rho_ref) = &findClustrMedoids_index();
	if($verbose){
		print  "classifycluster() ...\n";
	}
	my $cluster_ref = &Classifycluster($cluster_medoids_index_ref,$index_by_desc_rho_ref);
	if($verbose){
		print  "decide_clusterCore() ...\n";
	}
	my $halo_ref=&decide_clusterCore($cluster_ref,$rho_ref,$pdistances_ref);

	open(F,">$out_file");
	print F join("\n",@{$cluster_ref});
	close(F);

	open(F,">halo.txt");
	print F join("\n",@{$halo_ref});
	close(F);

	open(F,">cluster_medoids.txt");
	print F join("\n",@{$cluster_medoids_index_ref});
	close(F);
		

}

sub decide_clusterCore(){
	my ($cluster_ref,$rho_ref,$pdistances_ref) = @_;
	my @halo;
	my @border_rho;
	for my $index (0..(@{$cluster_ref}-1)){
		$border_rho[$index] = 0;
	}
	for my $key1(keys %{$pdistances_ref}){
		for my $key2 (keys %{$pdistances_ref->{$key1}}){
			if($pdistances_ref->{$key1}->{$key2} < $cutoff_distance && $cluster_ref->[$key1] != $cluster_ref->[$key2]){
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

sub findClustrMedoids_index(){
	if($verbose){
		print "step1.1 calculate_pdistance() ...\n";
	}
	my $pdistances_ref = &calculate_pdistance($cols_data_ref);
	if($verbose){
		print "step1.2 calculate_cutoff() ...\n";
	}
	$cutoff_distance = &calculate_cutoff($pdistances_ref,$t_value);
	if($verbose){
		print "cutoff_distance = $cutoff_distance...\n";
	}
	open(F,">cutoff_distance.txt");
	print F $cutoff_distance;
	close(F);
	if($verbose){
		print "step1.3 calculate_rho() ...\n";
	}
	my $rho_ref = &calculate_rho($pdistances_ref,$cutoff_distance);
	if($verbose){
		print "step1.4 get_index_by_desc_rho() ...\n";
	}
	my $index_by_desc_rho_ref = &get_index_by_desc_rho($rho_ref);
	my $deltas_ref = &calculate_delta($index_by_desc_rho_ref,$rho_ref);
	my $index_of_deltaRho_desc_ref=&calc_index_of_deltaRho_desc($rho_ref,$deltas_ref);
	my @cluster_medoids_index;
	if($k_value != -1){
		for(my $i = 0 ; $i < $k_value; $i++){
				$cluster_medoids_index[$i] = $index_of_deltaRho_desc_ref->[$i];
		}
	}else{
		my $cur_cluster = 1;
		$cluster_medoids_index[0] = $index_by_desc_rho_ref->[0];
		my $cur_deltaRho = $rho_ref->{$index_by_desc_rho_ref->[0]}*10 + $deltas_ref->{$index_by_desc_rho_ref->[0]};
		
		for(my $i = 1 ; $i < @{$index_by_desc_rho_ref} ; $i ++){
#			my $tmp_deltaRho = $rho_ref->{$index_by_desc_rho_ref->[$i]}*10 + $deltas_ref->{$index_by_desc_rho_ref->[$i]};
#			if($tmp_deltaRho > $cur_deltaRho / 2){
#				$cluster_medoids_index[$cur_cluster++] = $index_by_desc_rho_ref->[$i];
#			}
			if($rho_ref->{$index_by_desc_rho_ref->[$i]} > 0.6 && $deltas_ref->{$index_by_desc_rho_ref->[$i]} > 100){
				$cluster_medoids_index[$cur_cluster++] = $index_by_desc_rho_ref->[$i];
			}
		}
	}
	return (\@cluster_medoids_index,$rho_ref,$pdistances_ref,$index_by_desc_rho_ref);
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
#		$rho_ref->{$b} * ($deltas_ref->{$b}) <=> $rho_ref->{$a} * ($deltas_ref->{$a}) 
#		($rho_ref->{$b}) <=> ($rho_ref->{$a}) 
		($deltas_ref->{$b}) <=> ($deltas_ref->{$a}) 
	} keys %{$deltas_ref};
	return \@index_result;
}
sub calculate_delta{
	my ($index_by_rho_desc_ref,$rho_ref) = @_;
	my %deltas;
	my $max_dis = -1;
	for my $i (0..($data_num-1)){
		$deltas{$i} = 0;
	}
	
	for(my $i = 1 ; $i < @{$index_by_rho_desc_ref} ; $i ++){
		my $cur_delta = -1;
		my $cur_delta_end = -1;
		for(my $j = 0 ; $j < $i ; $j ++){
			my $distance = &claculate_ppdistance($index_by_rho_desc_ref->[$i],$index_by_rho_desc_ref->[$j]);
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
		print "delta_count:".$count."\n"
	}
	open(F,">rho_delta.txt");
	for my $key (sort {$a <=> $b} keys %{$rho_ref}){
		print F "$key\t".$rho_ref->{$key}."\t$deltas{$key}"."\n";
	}
	close(F);
	return \%deltas;
}
sub calculate_rho{
	my ($pdistances_ref,$cutoff_distance) = @_ ;
	my %rho_values;
	for (my $i=0 ; $i <$data_num ; $i++){
		$rho_values{$i} = 0;
	}
	foreach my $point_s (keys %{$pdistances_ref}){
		foreach my $point_e (keys %{$pdistances_ref->{$point_s}}){
			my $distance = $pdistances_ref->{$point_s}->{$point_e};
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
	my %pdistance;
	my $data_num = @{$cols_data_ref};
	for(my $i = 0 ; $i < $data_num-1 ; $i ++){
		for(my $j = $i+1 ; $j < $data_num; $j ++){
			my $dist = &claculate_ppdistance($i,$j);
			$pdistance{$i}{$j} = $dist;
		}
	}
	return \%pdistance;
}

sub claculate_ppdistance{
	my ($data_index_1,$data_index_2) = @_;
	my $line1_ref = $cols_data_ref->[$data_index_1];
	my $line2_ref = $cols_data_ref->[$data_index_2];
	my $distance = 0;
	for(my $i = 0 ; $i < @{$cols_ref} ; $i ++){
		if($line2_ref->[$i] eq ""){
			$line2_ref->[$i] = 0;
		}
		if($line1_ref->[$i] eq ""){
			$line1_ref->[$i] = 0;
		}
		$distance += ($line1_ref->[$i]-$line2_ref->[$i])**2;
	}
	$distance = sqrt($distance);
	return $distance;
}
sub calculate_cutoff{
	my ($pdistances_ref , $t_value) = @_;
	my @sorted_distance = &sort_distance($pdistances_ref);
	my $distance_num = scalar(@sorted_distance);
	if($distance_num != $data_num*($data_num-1)/2){
		print STDERR "data_num:".$data_num.",distance_num:" .$distance_num;
		print STDERR "distance not match data number!" and die;
	}
	my $k = &round($distance_num * $t_value);
	return $sorted_distance[$k-1];
}
sub sort_distance{
	my ($pdistances_ref) = @_ ;
	my @sorted_distance ;
	for my $k1 (keys %{$pdistances_ref}){
		for my $k2(keys %{$pdistances_ref->{$k1}}){
			$sorted_distance[@sorted_distance]= $pdistances_ref->{$k1}->{$k2};
		}
	}
	@sorted_distance = sort{$a <=> $b} @sorted_distance;
	return @sorted_distance;
}
