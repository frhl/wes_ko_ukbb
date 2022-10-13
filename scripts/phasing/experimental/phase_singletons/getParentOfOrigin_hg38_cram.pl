#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::HTS::Alignment;
use Bio::DB::HTS;
use Getopt::Long;

# This script is aim to detect the parent of origin for denovo (DN) mutations
# It checks the co-occurrence of DN mutation and the informative variants on reads or read pairs in the single sample
# Input text file contains variant pairs (DN mutation and infomrative variat pair)
# CAMs file is required
# Author: Shan Dong (2012)
# MOdified on 2017 -> worked on WGS
# Modified on 2019 -> swith to Bio::DB::HTS
# Modified on 2021 -> get options

my ($cram, $input_file, $input_prefix, $output_sum_file) =();
my $exit; 

GetOptions(
    'input|i=s' => \$input_file, # path for denovo indel list
    'cram|c=s' => \$cram, # path for cram files (separated by ";")
) or die "Usage: $0 --input (denovo mutation and informative variants pairs) --cram (cram file path) \n";

# print "-i:".$input_file."\n";
# print "-m:".$meta_file."\n";
# print "-b:".$cram_file."\n";

if ($input_file =~ m/^.*\/(.*)\.txt/){
    $input_prefix = $1;
    $output_sum_file = "Checking_results_summary.".$input_prefix.".txt";
}
else{
    print "Error: Wrong input variant file name format!\n";
}

################### Process the variant files######################
open (IN, "<".$input_file) or die "Can't open input file\n";
open (SUM, ">".$output_sum_file) or die "Can't open $!\n";

#ouput sum file header
#poo_result	confidence_level	distance	chr_dn	pos_dn	ref_dn	var_dn	sample	SSCID	fam	type	chr_info	pos_info	ref_info	var_info	info_origin	info_gt_code	dn_total_count	dn_support_count	percent_dn	info_total_count	info_support_count	percent_info	overlap_count	consist_count	percent_consist	incosist_count	percent_inconsist

print "processing:$input_file\n";
while (<IN>){
  chomp;
	# next if ($_ =~ m/^Chr/); #skip the header lines
	my $fix = $_;
    my @line = split/\s+/, $_;
    my $chr = $line[1];
    my $pos = $line[2];
    my $reference = $line[3];
    my $alt = $line[4];
    my $infopos = $line[10];
    my $inforef = $line[11];
    my $infovar = $line[12];
    my $info_poo = $line[13];

    if($info_poo eq "Ma"){
        $info_poo = "Mo"
    }

    my ($var_number,%var) = &get_info_from_bam($chr,$pos,$reference,$alt,$cram);
	my ($info_number,%info) =  &get_info_from_bam($chr,$infopos,$inforef,$infovar,$cram);

	###summary the results
	my $var_support_number =0;
	my $info_support_number =0;
	my $overlap_count =0;
	my $consist_count = 0;
	my $incosist_count =0;
	foreach my $key_info (sort keys %info){
		if ($info{$key_info} == 1){
		 $info_support_number = $info_support_number +1;
		}
	}
	foreach my $key_var (sort keys %var){
		if ($var{$key_var} == 1){
			$var_support_number = $var_support_number +1;
		}
		 foreach my $key_info (sort keys %info){
			 if($key_info eq $key_var){
				  $overlap_count = $overlap_count +1;
					if($var{$key_var} == 1 && $info{$key_info} == 1){
						$consist_count = $consist_count +1;
					}
					elsif($var{$key_var} == 0 && $info{$key_info} == 0){
						$consist_count = $consist_count +1;
					}
					elsif($var{$key_var} == 1 && $info{$key_info} == 0){
						$incosist_count = $incosist_count+1;
					}
					elsif($var{$key_var} == 0 && $info{$key_info} == 1){
						$incosist_count = $incosist_count+1;
					}
					else{
						next;
					}
			 }
	 	}
	}
	my ($percent_consist,$percent_inconsist);
	if($overlap_count == 0){
		$percent_consist =0;
		$percent_inconsist =0;
	}
	else{

		$percent_consist = sprintf ("%.2f",$consist_count/$overlap_count);
		$percent_inconsist = sprintf ("%.2f",$incosist_count/$overlap_count);
		if($consist_count == 0 ){
			$percent_consist =0;
		}
		if ($incosist_count == 0){
			$percent_inconsist =0;
		}
	}

  my ($percent_var,$percent_info);
	if ($var_number ==0 ){
			$percent_var =0;
	}
	elsif ($var_support_number == 0){
		 $percent_var =0;
	}
	else{
		$percent_var = sprintf ("%.2f",$var_support_number/$var_number);
	}
	if ($info_number ==0 ){
			$percent_info =0;
	}
	elsif ($info_support_number == 0){
		 $percent_info =0;
	}
	else{
		 $percent_info = sprintf ("%.2f",$info_support_number/$info_number);
	}

	## Start filtering by standard
	# Overlapped reads > 0
	# >10% of the reads support the DN variant
	# 10%-90% of the reads support the inforative variant
	# if Overlapped reads > 5, then envidence is strong; otherwise it's weak
	my $enviden;
	if ($overlap_count > 0){
		if ($percent_var > 0.1 && $percent_info > 0.1 && $percent_info <0.9){
			if($percent_consist > 0.9 && $percent_inconsist < 0.1){
				print SUM "$info_poo\t";
			}
			elsif($percent_inconsist > 0.9 && $percent_consist < 0.1){
				if ($info_poo eq 'Mo'){
					print SUM "Fa\t";
				}
				else{
					print SUM "Mo\t";
				}
			}
			else{
				print SUM "Unclear\t";
			}
			if ($overlap_count > 5){
				print SUM "Confident\t";
			}
			else{
				print SUM "Weak\t";
			}
			print SUM "$fix\t$var_number\t$var_support_number\t$percent_var\t";
			print SUM "$info_number\t$info_support_number\t$percent_info\t";
			print SUM "$overlap_count\t$consist_count\t$percent_consist\t$incosist_count\t$percent_inconsist\n";
		}
	}
}

close (IN);
close (SUM);

##########################################################################################
## Get information from the bam file using Bio::BD::HTS
## Give two variants, judge if they are exist in the same read or not
sub get_info_from_bam {
    my ($oldbamfile,$hts,$ref_len,$var_len,$thestart,$thestop,@alignments,@total,@thebases,$thebase,$read,$startpos,$sequence,$offset,$currentpos,$i,$maxlength,$maxseq,$cigar,$indelseq,$ref,$matched,$query,%gap,$sign,$thechange,$endpos,$ref_alignment,$var_alignment,$ref_count,$var_count,$t,$minlength,$matches,$analset,$real_ref_aln,$real_var_aln,$real_var_string,%count,%read_seq, $complex_reads,$rv_thechange,$rv_len);
	my %reads_varpos; #store the reads name which cover the variants, value =1 means contain the variant, otherwise no variant
    my ($thechromosome,$theposition,$thereference,$thevariant,$cramfile) =();
    ($thechromosome,$theposition,$thereference,$thevariant,$cramfile)= @_;

    if (-e $cramfile){
		$hts = Bio::DB::HTS->new(
        -bam => "$cramfile",
        #-autoindex
        );

        #$theposition is the start pos of an indel
        $thestart = $theposition;
        $ref_len = rindex ($thereference."\$","\$");
        $var_len = rindex ($thevariant."\$","\$");
        if ($ref_len == $var_len){
          $ref_alignment = $thereference;
          $var_alignment = $thevariant;
          $ref_alignment =~ tr/a-z/A-Z/;
          $var_alignment =~ tr/a-z/A-Z/;
					$thestop = $thestart;
					$sign = "*";
        }
        else{
          my $indel_size = abs ($ref_len - $var_len);
          if ($ref_len > $var_len) {
              #deletion
              $thestop = $thestart + $ref_len -1;
              $maxlength = $ref_len;
              $minlength = $var_len;
              $maxseq = $thereference;
              $sign = "-";
              $ref_alignment = $thereference; #for example: $ref_alignment is "CAAT"
              $var_alignment = $thevariant; #for example: $var_alignment is "C", so there is a deletion -ATT.
              for $i (1..($ref_len - $var_len)){
                  $var_alignment = $var_alignment."-"; #for example: $var_alignment is "C---" after finish this loop
              }
              $ref_alignment =~ tr/a-z/A-Z/;
              $var_alignment =~ tr/a-z/A-Z/;
          }
          else {
              #insertion
              $thestop = $thestart + $var_len -1;
              $maxlength = $var_len;
              $minlength = $ref_len;
              $maxseq = $thevariant;
              $sign = "+";
              $ref_alignment = $thereference;
              $var_alignment = $thevariant;
              for $i (1..($var_len - $ref_len)){
                  $ref_alignment = $ref_alignment."-";
              }
              $ref_alignment =~ tr/a-z/A-Z/;
              $var_alignment =~ tr/a-z/A-Z/;
          }

          $thechange = substr ($maxseq,$minlength,$maxlength-$minlength); #for example: $the_change here is "ATT"
          $thechange = $sign.$thechange; #for example: $the_change becomes "-ATT"
          $thechange =~ tr/a-z/A-Z/;
        }

        ###############################################################################################################

        #get alignment information: reads overlap with the specific indel here
        @alignments = $hts->get_features_by_location(
        -seq_id => $thechromosome,
        -start  => $thestart,
        -end    => $thestop,
        );

        if (@alignments){
            $complex_reads = 0; #record that suggest the position is a complex position, which given deletion event but see a insertion event
			my $read_count =0;
            #print "alignment\n";
            foreach $read (@alignments){
                $startpos = $read->start;     
                $endpos = $read->end;         
                my $q_start = $read->query->start;   
                my $q_end = $read->query->end;       
                my $read_name = $read->name;
                my $match_qual= $read->qual;
                if ($match_qual == 0){
                  next;
                }

                ($ref,$matches,$query)=$read->padded_alignment; #extract the alignmnet information
                my $len_r = rindex ($ref."\$","\$");
                my $len_m = rindex ($matches."\$","\$");
                my $len_q = rindex ($query."\$","\$");

                $cigar = $read->cigar_str;
                if ($cigar =~m/^(\d+)S/){
                    $ref = substr ($ref,$1,$len_r-$1);
                    $matches = substr ($matches,$1,$len_m-$1);
                    $query = substr ($query,$1,$len_q-$1);
                }
                #get the base count for each position in the indel region
                $offset = 0; #distance between indel position and read start
                $analset = 0; #distance between indel position and the read end
                $t = 0;
                if ($startpos < $thestart && $endpos > $thestop){
									  $read_count = $read_count +1;
										$reads_varpos{$read_name} = 0;
                    $offset = $thestart - $startpos; 
                    for my $i (0..$offset-1){
                        if (substr($ref,$i,1) eq "-"){ 
                            $t=$t+1;
                        }
                    }
                    $analset = $endpos -$thestart;
                    
                    #for SNP, sign = *
                    if ($sign eq "*"){
                            if (substr($query,$offset+$t,1) eq $var_alignment){
                                    $reads_varpos{$read_name} = 1;
                            }
                            elsif (substr($query,$offset+$t,1) eq $ref_alignment){
                                    $reads_varpos{$read_name} = 0;
                            }
                            else{
                                    $reads_varpos{$read_name} = 9;
                            }
                    }

                    #for deletion (sign eq "-") 
                    if ($sign eq "-"){
                        if (substr($query,$offset+$t,$maxlength) eq $var_alignment){ ##for deteltion exension or no change
                            $real_var_aln = substr($query,$offset+$t,1);
                            for my $x (1..$analset+1){
                                if (substr($query,$offset+$t+$x,1) eq "-"){
                                    $real_var_aln = $real_var_aln."-";
                                }
                                else{
                                    last;
                                }
                            }
                            if ($real_var_aln eq $var_alignment){
                                $reads_varpos{$read_name} = 1;
                            }
                            else{
                                $reads_varpos{$read_name} = 9;
                            }
                        }
                        elsif (substr($query,$offset+$t,$maxlength) =~ m/.*-.*[ATCGatcgnN]$/) { ##for deteltion shrink
                            $real_var_aln = substr($query,$offset+$t,1);
                            #print "first_base:$real_var_aln\n";
                            for my $x (1..$analset+1){
                                if (substr($query,$offset+$t+$x,1) eq "-"){
                                    $real_var_aln = $real_var_aln."-";
                                }
                                else{
                                    last;
                                }
                            }
                            $reads_varpos{$read_name} = 9;
                        }
                        elsif (substr($ref,$offset+$t,$maxlength) =~m/.*-.*/){
                            $reads_varpos{$read_name} = 9;
                        }
                        else{
                            $reads_varpos{$read_name} = 0;
                        }
                    }

                    #for insertion (sign eq "+") 
                    if ($sign eq "+"){
                        if(substr($ref,$offset+$t,$maxlength) eq $ref_alignment){ #for insertion exension or no change
                            $real_ref_aln = substr($ref,$offset+$t,1);
                            for my $x (1..$analset+1){
                                if (substr($ref,$offset+$t+$x,1) eq "-"){
                                    $real_ref_aln = $real_ref_aln."-";
                                }
                                else{
                                    last;
                                }
                            }
                            if ($real_ref_aln eq $ref_alignment){
                                $reads_varpos{$read_name} = 1;
                            }
                            else{
                                $reads_varpos{$read_name} = 9;
                            }
                        }
                        elsif (substr($ref,$offset+$t,$maxlength) =~ m/.*-.*/){ #for insertion shrink
                            $real_ref_aln = substr($ref,$offset+$t,1);
                            for my $x (1..$analset+1){
                                if (substr($ref,$offset+$t+$x,1) eq "-"){
                                    $real_ref_aln = $real_ref_aln."-";
                                }
                                else{
                                    last;
                                }
                            }
                            $reads_varpos{$read_name} = 9;
                        }
                        elsif (substr($query,$offset+$t,$maxlength) =~m/.*-.*/){
                            $reads_varpos{$read_name} = 9;
                        }
                        else{
                            $reads_varpos{$read_name} = 0;
                        }
                    }
                }
                else {
                    next;
                }
            }
			return ($read_count,%reads_varpos);
        }
    }
}
