#!/usr/bin/perl -w
use strict;
use Statistics::R;
my $R = Statistics::R->new();
$R->startR;
my @datafiles = (
	 "Data/flow_value_99_100.txt"
	,"Data/flow_value_114_115.txt"
	,"Data/flow_value2_99_100.txt"
	,"Data/flow_value3_99_100.txt"
	,"Data/flow_value3_103_104.txt"
);
-f "test_out" or mkdir("test_out",0755) ;

for (my $delta_min = 100; $delta_min >= 10 ; $delta_min-=10){
	my $data_count = @datafiles;
	$R->send(qq`jpeg(\"$delta_min\_0.02.jpg\",height=200*$data_count)`);
	$R->send(qq`par(mfrow=c($data_count,2));`);
	foreach my $data_file (@datafiles){
		if($data_file =~ /(flow_value\d?_(\d+)_(\d+))/){
			my $out_dir = "test_out\\$1";
			system("mkdir $out_dir");
			my $a = $2; my $b = $3;
			my $a_=$a-1; my $b_ = $b-1;
			for(my $i = 0.01 ; $i <=0.02 ; $i += 0.01){
				system("\"../density_cluster.pl\"  -file \"$data_file\" -o cluster.out -t_value $i -k 100 -delta_min $delta_min -rho_min 0 -g -cols $a_,$b_");

			system("Rscript testplot.R \"$data_file\" 0 $a $b");
			system("copy plot.jpg $out_dir\\$i.jpg");
			system("copy rho_delta.jpg $out_dir\\$i\_decide.jpg");
			system("copy cluster.out.cutoff_distance.txt $out_dir\\$i\_cutDis.txt");
			system("copy cluster.out.cluster $out_dir\\$i\_cluster.txt");
			system("copy cluster.out.cluster_medoids.txt $out_dir\\$i\_mediods.txt");
			if($i == 0.02){
					$R->send(qq`source('plotFunction.R') ; plotCluster(\"cluster.out\",\"$data_file\",0,$a,$b)`)
				}
			}

		}
	}
	$R->send(qq`dev.off();`);
	$R->stop();
}

