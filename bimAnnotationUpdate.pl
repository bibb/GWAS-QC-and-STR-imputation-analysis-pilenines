#!/usr/bin/perl
# updateLargeBIMAnnotation version 0.2
# compare two BIM files and get a new BIM with common SNPs from reference BIM
# BIM reference is Phase3 with multiallele SNPs
# 0.1: load firstly the query = SNPChip annotation 
#  and then check each element from reference sequently
# 0.2: double position error save in new format to stderr
#      bug: -9 korrekt, bei mehreren - ein auswaelen
# 0.3: only 1-26 as chromosome names
use strict;
use Data::Dumper;
# load query BIM
my %ids;
my %pos;
my @snps;
open(FIN, "< $ARGV[0]");
while(<FIN>) {
  chomp();
  my @a = split("\t");
  push @snps, \@a;
  # ids
  my $id_key = $a[1];
  print STDERR "ID_error\t$id_key\t$a[0]\t$a[3]\n" if(exists($ids{$id_key}));
  $ids{$id_key} = [];
  # pos
  my $chr = $a[0];
  $chr = 23 if($chr =~ /^X/);
  next if(!($a[0]>0) || $a[3]==0);
  my $pos_key = "$a[0]:$a[3]";
  print STDERR "position_error\t$id_key\t$a[0]\t$a[3]\n"
    if(exists($pos{$pos_key}));
  $pos{$pos_key} = [];
}
close(FIN);
# statistic
print STDERR "load chip SNPs:\t",scalar(@snps),"\n";
# load target BIM
while(<STDIN>) {
  chomp();
  my @a = split("\t");
  my $id_key = $a[1];
  if($id_key=~/^rs\d+:/){
    my @temp = split(":",$id_key);
    $id_key = $temp[0];
  }
  #^ bug phase 3 - biallele SNPs bestehen aus mehereren Basen
  if(length($a[4])>1 && length($a[4])==length($a[5])) {
    my $alleleA = substr($a[4],0,1);
    my $alleleB = substr($a[5],0,1);
    my $tempA = substr($a[4],1);
    my $tempB = substr($a[5],1);
    if($tempA eq $tempB) {
      $a[4] = $alleleA;
      $a[5] = $alleleB;
    }
  }
  #$ bug phase 3
  if(exists($ids{$id_key})) {
    push @{$ids{$id_key}}, \@a;
  }
  #^ chrX bug Mar2012
  my $chr = $a[0];
  $chr = 23 if($chr =~ /^X/);
  #$ chrX bug Mar2012
  my $pos_key = "$chr:$a[3]";
  if(exists($pos{$pos_key})) {
    push @{$pos{$pos_key}}, \@a;
  }
}
#print Dumper(\%ids),"\n";
#print Dumper(\%pos),"\n";

# compare annotations
foreach my $snp (@snps) {
  print join("\t",@$snp);
  # check by ID
  my $id_key = $$snp[1];
  my $ref_snps = $ids{$id_key};
  if(exists($ids{$id_key}) && scalar(@$ref_snps)>0) {
    print "\tid\t",choiseSNP($ref_snps, $snp),"\n";
    next;
  } 
  my $pos_key = "$$snp[0]:$$snp[3]";
  my $ref_snps = $pos{$pos_key};
  if(!exists($pos{$pos_key}) || scalar(@$ref_snps)==0) {
    print "\tno\t0\t0\n";
  } else { 
    print "\tpos\t",choiseSNP($ref_snps, $snp),"\n";
  }
}
#

# select correct SNP from list by chr,pos,allele1,allele2
sub choiseSNP {
  my $a = shift;
  my $query = shift;
  if(scalar(@$a)==1) {
    my $target = $$a[0];
    my $strandTest = check_strand($query, $target);
    return "$strandTest\t1\t".join("\t",@$target);
  }
  my $out = undef;
  my $outStrand = 0;
  my @b;
  foreach my $target (@$a) {
    next if($$query[0] ne $$target[0]); # chromosome
    next if($$query[3] ne $$target[3]); # position
    push @b, $target;
    my $strandTest = check_strand($query, $target);
    if($strandTest>0) {
      # next if($$query[0] !~ /^rs/); # wieso?
      die "error: duplicate SNP\t","\t",join(":",@$query),
        join(":",@$target),"\t",join(":",@$out) if($out != undef);
      $out = $target;
      $outStrand = $strandTest;
    }
  }
  if($outStrand>0) {
    return "$outStrand\t".scalar(@$a)."\t".join("\t",@$out);
  } 
  my $outText = "-9\t".scalar(@b);
  foreach my $target (@b) {
    $outText .= "\t".join("\t",@$target);
  }
  return $outText;
}

# check forward/reverse strand only for biallee.
sub check_strand {
 my $q = shift;
 my $r = shift;
 my $qA = $$q[4];
 my $qB = $$q[5];
 my $rA = $$r[4];
 my $rB = $$r[5];
 return 1 if($qA eq $rA && $qB eq $rB);
 return 2 if($qA eq $rB && $qB eq $rA);
# return 5 if($qA eq "0" && ($qB eq $rA || $qB eq $rB));
 $qA =~ tr/[ACGT]/[TGCA]/;
 $qB =~ tr/[ACGT]/[TGCA]/;
 return 3 if($qA eq $rA && $qB eq $rB);
 return 4 if($qA eq $rB && $qB eq $rA);
# return 6 if($qA eq "0" && ($qB eq $rA || $qB eq $rB));
 return -1;
}