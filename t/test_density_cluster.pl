#!/usr/bin/perl -w
use strict;
my $data_file = "Data/flow_value_122_125.txt";
mkdir("test_out",0755);
for(my $i = 0.01 ; $i <=0.1 ; $i += 0.01){
	system("\"../density_cluster.pl\"  -file \"$data_file\" -o cluster.out -t_value $i -g ");
	system("Rscript testplot.R $data_file 0");
#	system("copy Rplots.pdf test_out\\$i.pdf");
	system("copy plot.jpg test_out\\$i.jpg");
	system("copy rho_delta.jpg test_out\\$i\_decide.jpg");
	system("copy cutoff_distance.txt test_out\\$i\_cutDis.txt");
	system("copy cluster.out test_out\\$i\_cluster.txt");
	system("copy cluster_medoids.txt test_out\\$i\_mediods.txt");
	die;
}
