#! /usr/bin/perl

use strict;
use warnings;
use HTML::TableExtract;
use HTML::Table;
use Getopt::Long;
use File::Basename;

my ($runid, $html_out, $csv_out);
GetOptions("runid|R=s" =>\$runid,
	   "html|H=s"  =>\$html_out,
           "csv|C=s"   =>\$csv_out);

my $usage = qq{
USAGE: demux_summary.pl --runid /full/path/to/lane.html --html /full/path/to/out.html --csv /full/path/to/out.csv

Options :
	--runid	Full path to lane.html
        --html  Full path to output.html
	--csv   Full path to output.csv
};

die "$usage\n" unless ($runid && $html_out && $csv_out);

#my @runpath = split /\//, $runid;
#my @fields = split /_/, $runpath[$#runpath];
#my $flowcell;
#my $index = 0;
#foreach my $f (@fields) {
#	$index++;
#	last if ($f =~ m/\d\d\d\d/);
#} 
#$flowcell[$index] =~ s/^A//;
#my $path = "$runid/Reports/html/$flowcell[$index]/all/all/all";
my $path = dirname($runid);

open (OUT, ">$html_out") or die "Cannot open $html_out:$! \n";
open (OUT1, ">$csv_out") or die "Cannot open $csv_out:$! \n";

my %html = ( "Flowcell Summary" => "$path/lane.html", 
	     "Lane Summary"     => "$path/lane.html",
	     "Sample Summary"   => "$path/laneBarcode.html",
	     "Breakdown by Lane and Barcode" => "$path/laneBarcode.html");

my $out = [ "Flowcell Summary", "Lane Summary", "Sample Summary", "Breakdown by Lane and Barcode" ]; 
my %headers = ( "Flowcell Summary" => [ "Clusters \\(Raw\\)", "Clusters\\(PF\\)", "Yield \\(MBases\\)" ],
		"Lane Summary"     => [ 'Lane', 'PF Clusters' ],
		"Sample Summary"   => [ 'Lane', 'Project', 'Sample', 'Barcode sequence', 'PF Clusters' ],
		"Breakdown by Lane and Barcode" => [ 'Lane', 'Project', 'Sample', 'Barcode sequence', 'PF Clusters' ]
);
my $table_extract = HTML::TableExtract->new( depth => 0, count => 0); 
$table_extract->parse_file($html{"Flowcell Summary"});
print OUT1 "Flowcell ID: ";
foreach my $ts ($table_extract->tables) {
	my @row = $ts->rows;
	my @temp = split '\n', $row[0]->[0];
	$temp[0] =~ s/ \///;
	print OUT1 "\"", $temp[0], "\"\n\n";
	print OUT "<h1>Flowcell ID: $temp[0]</h1>\n";
}

foreach my $segment (@$out) {
	my $html = get_table($segment, $headers{$segment}, $html{$segment});
	print OUT1 "\n";
	if ($segment eq "Sample Summary") {
		print OUT "<h2>$segment</h2>\n";
		print OUT "<h5>(Per lane values averaged for - % of the lane, % Perfect barcode, % One mismatch barcode, % PF Clusters, % >= Q30 bases, Mean Quality Score)</h5>\n";
	}
	else { print OUT "<h2>$segment</h2>\n"; }
	print OUT "$html";
}

#system("mutt -e \"my_hdr Content-Type: text/html\" -s \"$runid demux summary\" `whoami`\@mail.nih.gov < $html_out");

sub get_table {

	my ($segment, $headers, $html) = @_;
	my $table_extract = HTML::TableExtract->new( headers => $headers, keep_headers => 1, slice_columns =>0 );
	my $tableOut = new HTML::Table();
#	$tableOut->setCaption($segment,"top");
	$tableOut->setBorder(1);
	$tableOut->setBGColor('papayawhip');
	$table_extract->parse_file($html);
	if ($segment eq "Sample Summary") {
		my %summary = ();
		my %count = (); #track numer of rows corresponding to each sample to calculate means.
		print OUT1 "$segment (Per lane values averaged for - % of the lane, % Perfect barcode, % One mismatch barcode, % PF Clusters, % >= Q30 bases, Mean Quality Score)\n";
		foreach my $ts ($table_extract->tables) {
			foreach my $row ($ts->rows) {
				foreach my $temp (@$row) {
					$temp =~ s/\n/ /; #remove line breaks in the headers.
					$temp =~ s/\,//g; 
				}
				if ($row->[0] eq "Lane") {
					shift @$row;
					shift @$row;
					print OUT1 "\"", join('","', @$row), "\"\n";
					$tableOut->addRow(@$row);
					$tableOut->setRCellsHead(1);
				}
				else { 
					if (exists $summary{$row->[2]}) {
						for(my $i=4; $i <= 11; $i++) {
							$summary{$row->[2]}[$i] += $row->[$i];
						}
						$count{$row->[2]} += 1;
					}
					else {
						$summary{$row->[2]} = [ @$row ];
						$count{$row->[2]} = 1;
					} 
				}
			}
		}
		foreach my $key (keys %summary) {
			my @index = (5, 6, 7, 9, 10 ,11);
			foreach my $i (@index) { 
				$summary{$key}[$i] = $summary{$key}[$i]/$count{$key};
			}
			shift @{ $summary{$key} };
			shift @{ $summary{$key} };	
			print OUT1 "\"", join('","', @{ $summary{$key} }), "\"\n";
			$tableOut->addRow(@{ $summary{$key} });
		}
	}
	else {	
		foreach my $ts ($table_extract->tables) {
			print OUT1 "$segment\n";
			foreach my $row ($ts->rows) {
				foreach my $temp (@$row) {
					$temp =~ s/\n/ /;
				}
				print OUT1 "\"", join('","', @$row), "\"\n";
				$tableOut->addRow(@$row);
				$tableOut->setRCellsHead(1);
  			}
		}
	}
	return($tableOut);
}
