#!/usr/bin/perl

#########################################################################
#
#  TarHunter
#
#  Version 1.0
#
#  This program is designed for predicting conserved miRNA target and 
#  target mimics in plants.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or    
#  (at your option) any later version.
#                                                                
#  This program is distributed in the hope that it will be useful,      
#  but WITHOUT ANY WARRANTY; without even the implied warranty of       
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       
#  GNU General Public License for more details.
#  
#  Author: Xuan Ma <xuanma@genetics.ac.cn>
#
#########################################################################

use strict;
use warnings;
use 5.010;
use Getopt::Long qw/:config no_ignore_case/;
use Cwd 'abs_path';
use File::Temp qw/tempfile tempdir/;
use File::Basename;
use IO::File;
use List::Util qw/max min/;

my $TarHunter_version = "1.0";
my $fasta_version   = "fasta36";

my $query_mir_list_file = '';
my $species_list_file = '';
my $dbs_folder = '';
my $output_folder = '';

my $program_option = 1;           #default: FASTA
my $mispair_cutoff = 100;         #set loose cutoff
my $seed_mispair_cutoff = 100;    #set loose cutoff
my $mimics_arm_misp_cutoff = 0;   #default: no misp in eTM arms
my $bits_cutoff = 0.1;            #default: 10% of best bitscore
my $nc_len = 100;                 #default: homo_mode nc/CDS 100 bp
my $iden_cutoff = 0.7;            #default: homo_mode identity 70%
my $jobs = 1;                     #default: 1 job
my $threads = 1;                  #default: 1 CPU

my $nonmirbase_file;
my $MUSCLE_alignment_file;
my $score_cutoff;
my $no_orth_mir_search;
my $no_conservation;
my $nc_target;
my $mimics;
my $tab_format;
my $help;


GetOptions (
    "n|non_mirbase=s"  => \$nonmirbase_file,
    "q|qmir=s"         => \$query_mir_list_file,
    "s|spe=s"          => \$species_list_file,
    "b|db=s"           => \$dbs_folder,
    "p|prog=i"         => \$program_option,
    "a|aln=s"          => \$MUSCLE_alignment_file,
    "M|total_misp=i"   => \$mispair_cutoff,
    "m|seed_misp=i"    => \$seed_mispair_cutoff,
    "f|score=f"        => \$score_cutoff,
    "c|iden=f"         => \$iden_cutoff,
    "N|no_conserve"    => \$no_conservation,
    "R|no_ortho_mir"   => \$no_orth_mir_search,
    "G|nc_targ"        => \$nc_target,
    "I|mimics"         => \$mimics,
    "i|mimics_str=i"   => \$mimics_arm_misp_cutoff,
    "L|len=i"          => \$nc_len,
    "j|jobs=i"         => \$jobs,
    "T|threads=i"      => \$threads,
    "t|tab"            => \$tab_format,
    "o|output=s"       => \$output_folder,
    "h|help"           => \$help
);


$score_cutoff = 100 if $mimics;            #eTM: set loose score cutoff
$score_cutoff = 4 if (! $mimics and ! $score_cutoff); #default score: 4


STDERR->say ("\n", "=="x40);
STDERR->say ("TarHunter(version $TarHunter_version)    Plant miRNA target searching tool");
STDERR->say ("=="x40);

my $time = scalar localtime;
STDERR->say ("$time");

my @cmdline = ();
my $buffer;
if (open my $cmd_fh, "<:raw", "/proc/$$/cmdline") {
    read($cmd_fh, $buffer, 8388608); 
    close $cmd_fh;
    @cmdline = split /\0/s, $buffer;
}
STDERR->say ("Command Line: ", join(" ", @cmdline), "\n");

STDERR->say ("-"x6, " Checking arguments and tools ...");

#check arguments
check_arguments();


#check if output folder already exists
if (-e $output_folder) {
    STDERR->say ("Output folder already exists. TarHunter exits.");
    exit;
}


#get TarHunter directory path
my $program_path = abs_path($0);
my $dir = dirname($program_path);


#create temporary folder
system "mkdir $dir/temp" unless -e "$dir/temp";
my $temp_dir = tempdir('TarHunter_job_XXXXXXXXX', DIR => "$dir/temp");
STDERR->say ("Creating temporary folder ($temp_dir).");

system "mkdir $output_folder $temp_dir/cdhit $temp_dir/ortho $temp_dir/cdsnew $temp_dir/pr_db $temp_dir/mimics";

my $usearch_version = "$dir/bin/usearch8.1.1756_i86linux32";
my $muscle_version  = "$dir/bin/muscle3.8.31_i86linux64";
system "chmod a+x $usearch_version" unless -x $usearch_version;
system "chmod a+x $muscle_version" unless -x $muscle_version;

#check tools USEARCH, MUSCLE, MCL, FASTA and RNAhybrid
check_tools();



####################  Reading miRNAs  ####################

#read mir list and species list
STDERR->say ("\n", "-"x6, " Reading miRNA file ...");

if ($nonmirbase_file) {
    system "cat $nonmirbase_file $dir/miRs/miRBase.fa > $temp_dir/cdhit/miRBase.fa";
}

my $mir_list_href = read_list($query_mir_list_file);
my $species_list_href = read_list($species_list_file);
my $num_of_species = scalar (keys %{ $species_list_href});


#produce new miR list and new miR cluster file
if ($nonmirbase_file) {
    system "cd-hit-est -i $temp_dir/cdhit/miRBase.fa -o $temp_dir/cdhit/out -c 0.85 -n 5 -d 40 -g 1 >>$temp_dir/msg.txt 2>&1";
    cdhit_parser();
}

my %mir_newlist = %{$mir_list_href};
orth_mir() unless $num_of_species == 1;


#read mirBase sequence file
my $miRBase;
if ($nonmirbase_file) {$miRBase = "$temp_dir/cdhit/miRBase.fa";}
else                  {$miRBase = "$dir/miRs/miRBase.fa";}
my $mir_seq_href = read_seq_file($miRBase);


#get sequences of new miR list
my %mir_newseq;
for my $id (keys %{$mir_seq_href}) {
    if (exists $mir_newlist{$id}) {
        $mir_newseq{$id} = $mir_seq_href->{$id};
    }
}

if (! %mir_newseq) {
    STDERR->say ( "No miRBase miRNAs. TarHunter exits." );
    exit;
}
undef $mir_seq_href;


 
####################  Orthologous gene clustering  ####################

my ($MUSCLE_aln_href, $all_proteomes_href);
unless (defined $no_conservation or $num_of_species == 1 or
        defined $mimics or defined $nc_target or defined $MUSCLE_alignment_file) {
        
    STDERR->say ("\n", "-"x6, " Orthologous gene clustering ...");
    
    #convert cds dbs to protein dbs
    $all_proteomes_href = generate_pr_dbs();

    #run UBLAST, MCL, MUSCLE
    my $gene_group_href;
    $gene_group_href = run_UBLAST_MCL();
    run_MUSCLE($gene_group_href);
}



####################  miRNA target search  ####################

if ($program_option == 0) {
    STDERR->say ("\n", "-"x6, " TarHunter analysis completed.");
    exit;
}

STDERR->say ("\n", "-"x6, " miRNA target search ...");

#read muscle alignment file
if ( defined $MUSCLE_alignment_file ) {
    $MUSCLE_aln_href = read_seq_file($MUSCLE_alignment_file);
}


#run FASTA and RNAhybrid for each miR
my ($rnahy_output_ID, $fasta_output_ID) = (0, 0);

for my $query_mir (sort keys %mir_newseq) {
    next unless $query_mir =~ /^(...)\-/;
    my $cds_file;
    my @db_files = glob "$dbs_folder/*.fa";
    next unless ($cds_file) = grep { /$1\.fa/ } @db_files;
    $cds_file = basename $cds_file;

    my $mir_tmpfile = IO::File->new("$temp_dir/mir.fa.tmp", 'w') or die $!;
    $mir_tmpfile->say (">$query_mir\n$mir_newseq{$query_mir}");
    $mir_tmpfile->close;

    if ($program_option == 1) {
        run_FASTA($query_mir, "$temp_dir/mir.fa.tmp", "$dbs_folder/$cds_file");
        
    } elsif ($program_option == 2) {
        unless ( -e "$temp_dir/cdsnew/$cds_file\.cdsnew" ) {

            #read cds db
            my $cds_file_href =  read_seq_file("$dbs_folder/$cds_file");
                
            #creat new cds db for RNAhybrid
            STDERR->say ("Producing new $cds_file file for RNAhybrid...");
            my $cdsnew_href = cds_new($cds_file_href);

            my $cdsnew_file = IO::File->new("$temp_dir/cdsnew/$cds_file\.cdsnew", 'w') or die $!;
            foreach (keys %{$cdsnew_href}) {
                $cdsnew_file->say (">$_\n$cdsnew_href->{$_}");
            }
        }
        run_RNAhybrid($query_mir, "$temp_dir/mir.fa.tmp", "$temp_dir/cdsnew/$cds_file\.cdsnew");
    }
}



####################  Conservation filter  ####################

if ($no_conservation or $num_of_species == 1) {
    STDERR->say ("\n", "-"x6, " TarHunter analysis completed."); 
    exit;
}

if ($mimics or $nc_target) { # nc/mimics target conservation
    STDERR->say ("\n", "-"x6, " Finding cross-species conserved target sites ...");
    mimics_nctarg_conservation();
    
} else { # CDS target conservation
    STDERR->say ("\n", "-"x6, " Finding cross-species conserved CDS target sites ...");
    system "sort -k 2,2 $temp_dir/mircluster.tmp > $temp_dir/mircluster.sorted.tmp";

    my $program;
    if ($program_option == 1)    { $program = 'fasta'; }
    elsif ($program_option == 2) { $program = 'rnahy'; }

    conservation("$temp_dir/$program\_out", "$temp_dir/$program\_conserv_out");
    undef $MUSCLE_aln_href;

    if ($tab_format) {
        system "head -n 1 $temp_dir/$program\_conserv_out >> $output_folder/$program\_out_conserv";
        system "tail -n +2 $temp_dir/$program\_conserv_out | sort -u >> $output_folder/$program\_out_conserv";
    } else {
        remove_duplicate("$temp_dir/$program\_conserv_out", "$output_folder/$program\_out_conserv");
    }
}

system "rm -r $temp_dir/cdhit $temp_dir/ortho $temp_dir/cdsnew $temp_dir/pr_db $temp_dir/mimics";
system "rm $temp_dir/*.tmp";
system "rm $temp_dir/*_out" if -e "$temp_dir/*_out";

STDERR->say ("\n", "-"x6, " TarHunter analysis completed.\n");



####################  Subroutines  ####################

sub check_arguments {
    unless ($query_mir_list_file) {
        STDERR->say ("FATAL: query miR list file does not exist.");
        usage();
    }
    unless ($dbs_folder) {
        STDERR->say ("FATAL: dbs folder is not provided.");
        usage();
    }
    unless ($species_list_file) {
        STDERR->say ("FATAL: species list file does not exist.");
        usage();
    }
    unless ($output_folder) {
        STDERR->say ("FATAL: output folder is not provided.");
        usage();
    }
    if ($MUSCLE_alignment_file and !(-e $MUSCLE_alignment_file)) {
        STDERR->say ("FATAL: protein alignment file does not exists.");
        usage();
    }
    if ($nonmirbase_file and !(-e $nonmirbase_file)) {
        STDERR->say ("FATAL: non-miRBase miRNA file does not exists.");
        usage();
    }
    if ($help) { usage(); }
}


sub usage {
my $usage = << "USAGE";

Usage:  
    perl TarHunter.pl -q <mir_list> -s <spe_list> -b <dbs_dir> -o <output_dir> [Options]

Required arguments:
    -q (--qmir):         miRNA list
    -s (--spe):          species list
    -b (--db):           dbs directory
    -o (--output):       output directory

Options:
    -p (--prog):         programs
                         (1:FASTA, 2:RNAhybrid, 0:neither) [Default:  1 ]
    -a (--aln):          protein alignment file            [Default: off]
    -n (--nonmirbase):   non-miRBase miRNA file            [Default: off]  

    -M (--total_misp):   max. total mispairs               [Default: off]
    -m (--seed_misp):    max. seed mispairs                [Default: off]
    -f (--score):        score cutoff                      [Default:  4 ]

    -G (--nc_targ):      homo_mode conservation filter     [Default: off]
    -c (--iden):         homo_mode target site identity    [Default: 0.7]
    -L (--len):          homo_mode target site length      [Default: 100]

    -I (--mimics):       eTM search                        [Default: off]
    -i (--mimics_str):   eTM stringency
                         (0: strict, 1: relaxed)           [Default:  0 ]

    -N (--no_conserve):  no ortho_mode conservaton filter  [Default: off]
    -R (--no_ortho_mir): no orthologous miRNA search       [Default: off]

    -j (--jobs):         jobs                              [Default:  1 ]
    -T (--threads):      threads                           [Default:  1 ]
    -t (--tab):          tabular format output             [Default: off]
    -h (--help):         help information

USAGE

    print $usage;
    exit;
}



sub check_tools {
    my $CDHIT_command = qx(which cd-hit-est 2>>$temp_dir/msg.txt);
    unless($CDHIT_command) {
        STDERR->say ("Warning: cd-hit-est not found.");
    }
    
    my $MCL_command = qx(which mcl 2>>$temp_dir/msg.txt);
    unless($MCL_command ) {
        STDERR->say ("Warning: MCL not found.");
    }
    
    my $PARALLEL_command = qx(which parallel 2>>$temp_dir/msg.txt);
    unless($PARALLEL_command ) {
        STDERR->say ("Warning: GNU parallel not found.");
    }
    
    my $RNAhybrid_command = qx(which RNAhybrid 2>>$temp_dir/msg.txt);
    unless ($RNAhybrid_command) {
        STDERR->say ("Warning: RNAhybrid not found.");
    }
    
    my $FASTA_command = qx(which $fasta_version 2>>$temp_dir/msg.txt);
    unless ($FASTA_command) {
        STDERR->say ("Warning: $fasta_version not found.");
    }
}


sub read_list {
    my $file = shift;
    my %list;
    my $infile = IO::File->new( $file, 'r' ) or die $!;
    while (<$infile>) {
        chomp; 
        s/\r//;
        s/^\s*|\s*$//;
        next if /^\s*$/;
        $list{$_} = 1;
    }
    return \%list;
}


sub read_seq_file {
    my ($seqfile) = @_;
    my %seqhash;

    local $/ = '>';
    my $infile = IO::File->new($seqfile, 'r') or die $!;
    while (<$infile>) {
        s/\r//g;
        chomp;
        my ($id, $seq) = /(\S+)(?:.*?)\n(.+)/s or next;
        ($seq = uc $seq) =~ s/\n//g;
        $seqhash{$id} = $seq;
    }
    return \%seqhash;
}


sub orth_mir {
    my $mir_group_href;
    if ($nonmirbase_file) { $mir_group_href = ortho_group("$temp_dir/cdhit/mir_cluster.txt") }
    else                  { $mir_group_href = ortho_group("$dir/miRs/mir_cluster.txt") }
    
    #output mircluster.tmp file for conservation analysis
    my $mircluster_tmpfile = IO::File->new("$temp_dir/mircluster.tmp", 'w') or die $!;
    my %all_spe;
    MIR:for my $mir (keys %{$mir_list_href}) {
        my %passed_spe;
        next unless $mir =~ /^(...)-miR/;
        $passed_spe{$1} = 1;
        $all_spe{$1} = 1;
        
        #find group ID of query mir 
        for my $group (keys %{$mir_group_href}) {
            next unless exists $mir_group_href->{$group}{$mir};
            $mircluster_tmpfile->say ("$group\t$mir");
            
            #without orthologous mir search, only output mircluster.tmp file
            next MIR if defined $no_orth_mir_search;
            
            #find the most orthologous mir from the file mir_cluster.txt
            local $/ = '#';
            my $mircluster_file;
            if ($nonmirbase_file) { $mircluster_file = IO::File->new("$temp_dir/cdhit/mir_cluster.txt", 'r') or die $! }
            else                  { $mircluster_file = IO::File->new("$dir/miRs/mir_cluster.txt", 'r') or die $! }
    
            my $omitted = <$mircluster_file>;
            while (<$mircluster_file>) {
                chomp; 
                s/\r//g;
                next unless /^$group\n/;
                my @mir_arr = split/\n/, $_;
            
                for my $mir_member (@mir_arr) {
                    next unless $mir_member =~ /^(...)-miR/;
                    next if exists $passed_spe{$1};
                    next unless exists $species_list_href->{$1};
                    $mir_newlist{$mir_member} = 1;
                    $mircluster_tmpfile->say ("$group\t$mir_member");
                    $passed_spe{$1} = 1;
                }
            }
        }
    }
    if (! %all_spe) { STDERR->say("No miRBase miRNAs found. TarHunter exits."); exit; }
    $num_of_species == 1 if scalar keys %all_spe == 1;
}


sub generate_pr_dbs {
    my %all_proteomes;
    for my $species (keys %{$species_list_href}) {
        my $cds_file;
        my %proteome;
        my @db_files = glob "$dbs_folder/*.fa";
        next unless ($cds_file) = grep { /$species\.fa/ } @db_files;
    
        (my $pr_file = basename $cds_file) =~ s/\.fa/\.pr\.fa/;
    
        #read cds db
        my $cds_href = read_seq_file($cds_file);
    
        #translate cds db
        my $pr_outfile = IO::File->new("$temp_dir/pr_db/$pr_file", 'w') or die $!;
        for my $id (keys %{$cds_href}) {
            $proteome{$id} = translation($cds_href->{$id});
            $pr_outfile->say (">$id\n$proteome{$id}");
        }
        $pr_outfile->close;
    
        #merge all pr dbs
        STDERR->say ("Converting $species cds dbs to protein dbs ...");
        @all_proteomes{keys %proteome} = values %proteome;
    }
    return \%all_proteomes;
}


sub cdhit_parser {
    my ($k, $m, $v);
    my %h;

    open my $in_fh, '<', "$temp_dir/cdhit/out.clstr" or die $!;
    while(<$in_fh>){
        chomp;
        if (/>Cluster/) {
            s/>/#/; s/\s+//;
            $k = $_;
        } else {
            next unless />(...-miR.+)\.\.\..*(\/.+%|\*)$/;
            $m = $1;
            $v = $2 =~ s/%|\///gr;
            $v = 0 if $v eq '*';
            push @{ $h{$k} }, [$m, $v];
        }
    }
    close $in_fh;

    open my $out, '>', "$temp_dir/cdhit/mir_cluster.txt" or die $!;
    for my $i (sort keys %h) {
        print $out "$i\n";
        for my $j ( @{ $h{$i} } ) {
            next unless $j->[1] == 0;
            print $out "$j->[0]\n";
            last;
        }
        for my $l ( sort {$b->[1] <=> $a->[1]} @{ $h{$i} } ) {
            next if $l->[1] == 0;
            print $out "$l->[0]\n";
        }
    }
    close $out;
}


sub translation {
    my $DNA = shift;
    $DNA = uc $DNA;

    my %genetic_code = (
    'ATA' => 'I',    'ATC' => 'I',    'ATT' => 'I',    'ATG' => 'M',    'ACA' => 'T',
    'ACC' => 'T',    'ACG' => 'T',    'ACT' => 'T',    'AAC' => 'N',    'AAT' => 'N',
    'AAA' => 'K',    'AAG' => 'K',    'AGC' => 'S',    'AGT' => 'S',    'AGA' => 'R',
    'AGG' => 'R',    'TCA' => 'S',    'TCC' => 'S',    'TCG' => 'S',    'TCT' => 'S',
    'TTC' => 'F',    'TTT' => 'F',    'TTA' => 'L',    'TTG' => 'L',    'TAC' => 'Y',
    'TAT' => 'Y',    'TAA' => '*',    'TAG' => '*',    'TGC' => 'C',    'TGT' => 'C',
    'TGA' => '*',    'TGG' => 'W',    'CTA' => 'L',    'CTC' => 'L',    'CTG' => 'L',
    'CTT' => 'L',    'CCA' => 'P',    'CCC' => 'P',    'CCG' => 'P',    'CCT' => 'P',
    'CAC' => 'H',    'CAT' => 'H',    'CAA' => 'Q',    'CAG' => 'Q',    'CGA' => 'R',
    'CGC' => 'R',    'CGG' => 'R',    'CGT' => 'R',    'GTA' => 'V',    'GTC' => 'V',
    'GTG' => 'V',    'GTT' => 'V',    'GCA' => 'A',    'GCC' => 'A',    'GCG' => 'A',
    'GCT' => 'A',    'GAC' => 'D',    'GAT' => 'D',    'GAA' => 'E',    'GAG' => 'E',
    'GGA' => 'G',    'GGC' => 'G',    'GGG' => 'G',    'GGT' => 'G',
    );

    my $protein;
    my $i = 0;
    while ( $i <= length($DNA) - 3 ) {
        my $aa = substr($DNA, $i, 3);

        if ( exists $genetic_code{$aa} ) {
            $protein .= $genetic_code{$aa};
        } else {
            $protein .= 'X';
        }
        $i += 3;
    }
    $protein = $1 if $protein =~ /(.+?)\*/; # remove AAs behind 1st stop codon
    return $protein;
}


sub run_UBLAST_MCL {
    STDERR->say ("Running UBLAST ...");
    
    my $db_list = IO::File->new ("$temp_dir/pr_db/db_list.txt", 'w') or die $!;

    #make ublast db (reciprocal ublast method)
    my @pr_dbs = glob "$temp_dir/pr_db/*.pr.fa";
    for my $db1 (@pr_dbs) {
        (my $db1_name = basename $db1) =~ s/\.pr\.fa//;
        system "$usearch_version -makeudb_ublast $db1 -output $temp_dir/pr_db/$db1_name.udb >>$temp_dir/msg.txt 2>&1";
        
        #query each db
        for my $db2 (@pr_dbs) {
            next if $db1 eq $db2;
            (my $db2_name = basename $db2) =~ s/\.pr\.fa//;
            $db_list->print ("$db2\t$temp_dir/pr_db/$db1_name.udb\t$temp_dir/pr_db/$db2_name\_$db1_name\_result\n");
        }
    }
    $db_list->close;
    
    system "parallel --no-notice --jobs $jobs --colsep '\t' -a $temp_dir/pr_db/db_list.txt $usearch_version -ublast {1} -db {2} -threads $threads -evalue 1e-5 -top_hit_only -id 0.5 -blast6out {3} >>$temp_dir/msg.txt 2>&1";

    # reciprocal best hits (RBH) 
    my (%passed, %reciprocal_pairs, %gene_group);
    my @res_files = glob "$temp_dir/pr_db/*result";
    for my $file1 (@res_files) {
        if ($file1 =~ /pr\_db\/(.+?)\_(.+?)\_result/) {
            my ($spe1, $spe2) = ($1, $2);
            next if exists $passed{"$spe2 $spe1"};
            $passed{"$spe1 $spe2"} = 1;
            (my $file2) = grep { /pr_db\/$spe2\_$spe1\_result/ } @res_files or next;
            
            my $pair1 = protein_pairs($file1);
            my $pair2 = protein_pairs($file2);
            
            for (keys %{$pair1}) {
                next unless /(.+)\s(.+)/;
                if (exists $pair2->{"$2 $1"}) {
                    $reciprocal_pairs{"$_"} = 
                        $pair1->{$_} . " " . $pair2->{"$2 $1"}; # get reciprocal hits
                }
            }
        }
    }
    
    my (%bits, %best_bits);
    for (keys %reciprocal_pairs) {
        my ($gene1, $gene2) = split;
        my ($eval1, $bits1, $eval2, $bits2) = split /\s/, $reciprocal_pairs{$_};
        push @{$bits{$gene1}}, $bits1;
        push @{$bits{$gene2}}, $bits2;
    }
    for (keys %bits) { $best_bits{$_} = max @{$bits{$_}}; } # get best bitscore
    
    for (keys %reciprocal_pairs) {
        my ($gene1, $gene2) = split;
        my ($eval1, $bits1, $eval2, $bits2) = split /\s/, $reciprocal_pairs{$_};
        if ($bits1 / $best_bits{$gene1} < $bits_cutoff or
            $bits2 / $best_bits{$gene2} < $bits_cutoff) { # bitscore cutoff
            delete $reciprocal_pairs{$_}; 
        }
    }
    
    # MCL orthologous gene clustering
    STDERR->say ("Running MCL ...");
    my $abc = IO::File->new ("$temp_dir/pr_db/ublast.abc", 'w') or die $!;
    for (keys %reciprocal_pairs) {
        my ($gene1, $gene2) = split;
        my ($eval1, $bits1, $eval2, $bits2) = split /\s/, $reciprocal_pairs{$_};
        print $abc "$gene1\t$gene2\t$bits1\n$gene2\t$gene1\t$bits2\n";
    }
    $abc->close;
    
    system "mcxload -abc $temp_dir/pr_db/ublast.abc --stream-mirror -o $temp_dir/pr_db/ublast.mci -write-tab $temp_dir/pr_db/ublast.tab 2>>$temp_dir/msg.txt";
    system "mcl $temp_dir/pr_db/ublast.mci -I 1.4 -use-tab $temp_dir/pr_db/ublast.tab -o $temp_dir/pr_db/ublast.out 2>>$temp_dir/msg.txt";

    my $group_num;
    my $mcl_result = IO::File->new ("$temp_dir/pr_db/ublast.out", 'r') or die $!;
    while (<$mcl_result>) {
        chomp;
        next if /^\s*$/;
        my $group = 'Group' . ++$group_num;
        my @gene_arr = split;
        $gene_group{$group}{$_} = 1 for (@gene_arr);
    }
    return \%gene_group;
}


sub protein_pairs {
    my $infile = shift;
    my %result;
    my $infile_fh = IO::File->new ($infile, 'r') or die $!; 
    while (<$infile_fh>) {
        my @arr = split;
        
        if ($result{"$arr[0] $arr[1]"}) { # query, subj
            my ($new_evalue, $new_bitscore) = split /\s/, $result{"$arr[0] $arr[1]"};
            next if ($arr[-2] > $new_evalue or $arr[-1] < $new_bitscore);
        }
        $result{"$arr[0] $arr[1]"} = "$arr[-2] $arr[-1]";
    }
    return \%result;
}


sub ortho_group {
    my $file = shift;
    my %ortho;
    my $group;

    my $infile = IO::File->new( $file, 'r' ) or die $!;
    while(<$infile>) {
        next if /^\s+$/;
        chomp;
        s/\r//;
        if (/^#(\w+)/) {
            $group = $1;
            next;
        }
        $ortho{$group}{$_} ++;
    }
    return \%ortho;
}


sub run_MUSCLE {
    my $gene_group_href = shift;
    
    #generate individual orthologous groups;
    STDERR->say ("Generating othologous protein groups ...");

    for my $groupID (keys %{$gene_group_href}) {
        for my $member ( keys %{$gene_group_href->{$groupID}} ) {
            next unless exists $all_proteomes_href->{$member};
            if ($all_proteomes_href->{$member} =~ /\*/) { #omit prs with stop-codons
                STDERR->say ("$member translated pr (+1 frame) contains stop codon, omitted.");
                next;
            }
            
            open my $group_file, '>>', "$temp_dir/ortho/$groupID.fa";
            print $group_file ">$groupID\|$member\n", $all_proteomes_href->{$member}, "\n";
            close $group_file;
        }
    }
    
    undef $all_proteomes_href;
    
    #running MUSCLE
    STDERR->say ("Running MUSCLE ...");
    
    system "find $temp_dir/ortho/ -name '*.fa' | parallel --no-notice --jobs $jobs $muscle_version -quiet -in {} -out {.}.afa";
    
    # store MUSCLE results into hash and output to ortho_aln file
    my $AFA_file = IO::File->new("$output_folder/ortho_aln_MUSCLE.afa", 'w') or die $!;
    my $MUSCLE_file = IO::File->new("find $temp_dir/ortho/ -name '*.afa' | xargs cat |") or die $!;
    my $id  = '';
    my $seq = '';
    while (<$MUSCLE_file>) {
        chomp;
        if (/^>(\S+)/) {
            if ($id ne '') {
                $MUSCLE_aln_href->{$id} = $seq;
                $AFA_file->say (">$id\n$seq");
            }
            $id  = $1;
            $seq = '';
        } else {
            $seq .= $_;
        }
    }
    $MUSCLE_aln_href->{$id} = $seq;
    $AFA_file->say (">$id\n$seq");
    
    STDERR->say ("MUSCLE completed.");
}


sub cds_new { #cds are split into 2kb segments, as RNAhyb corrupts with long seqs
    my $cds_file_href = shift;
    my %cdsnew;

    while ( my ($id, $seq) = each %{$cds_file_href} ) {

        my $line_num = 2 * int(length($seq)/2000);

        my $new_line_num;
        if (length($seq) <= 2000) {
            $new_line_num = 1;
        } elsif (length($seq)%2000 < 15) {
            $new_line_num = $line_num;
        } else {
            $new_line_num = $line_num + 1;
        }

        for my $i (1..$new_line_num) {
            if ($i%2) {
                my $start  = 1000 * ($i - $i%2) + 1;
                my $newseq = substr($seq, $start-1, 2000);
                $cdsnew{"$id\_$start"} = $newseq;

            } else {
                my $start  = 1000 * ($i - $i%2) - 14;
                my $newseq = substr($seq, $start-1, 30);
                $cdsnew{"$id\_$start"} = $newseq;
            }
        }
    }
    return \%cdsnew;
}


sub run_RNAhybrid {
    my ($mirID, $mirfile, $cdsnewfile) = @_;
    
    STDERR->print ("Running RNAhybrid for $mirID ...  ");

    open my $run_rnahyb, "RNAhybrid -d theta -q $mirfile -t $cdsnewfile -u $mispair_cutoff -v $mispair_cutoff -c |" or die $!;
    while (<$run_rnahyb>) {
        $rnahy_output_ID ++;
        chomp;
        RNAhybrid_parser($_);
    }
    STDERR->say ("done.");
}


sub RNAhybrid_parser {
    my $line = shift;
    my @line_arr = split/:/,$line or return -1;
    my ($field7, $field8, $field9, $field10) = @line_arr[7, 8, 9, 10];

    return -1 unless ($field7 and $field8 and $field9 and $field10);
    my $field7_copy  = $field7  =~ s/\s+//gr;
    my $field10_copy = $field10 =~ s/\s+//gr;
    if (length($field7_copy) - 2 > $mispair_cutoff or length($field10_copy) > $mispair_cutoff) {
        return -1;
    }

    my @f7  = split//, reverse $field7;
    my @f8  = split//, reverse $field8;
    my @f9  = split//, reverse $field9;
    my @f10 = split//, reverse $field10;

    #find mir 5' end
    my ($mir_beg, $mir_end);
    for my $i (0..$#f7) {
        next if ($f9[$i]=~/\s/ and $f10[$i]=~/\s/);
        $mir_beg = $i;
        last;
    }
    
    #find mir 3' end
    for my $i (-$#f7..0) {
        next if ($f9[-$i]=~/\s/ and $f10[-$i]=~/\s/);
        $mir_end = -$i;
        last;
    }

    my $cleavage = 'Yes';
    my ($tar, $mir, $aln_symb);
    my ($pos_count, $score, $g3_g5, $mispair, $seed_mispair, $p8_mispair,
        $middle_mispair, $end_mispair, $MFE_ratio) = (0,0,0,0,0,0,0,0,0,0);

    for my $i ($mir_beg..$mir_end) {
        my $tar_loop_nt = $f7[$i];
        my $tar_pair_nt = $f8[$i];
        my $mir_pair_nt = $f9[$i];
        my $mir_loop_nt = $f10[$i];

        if ($tar_loop_nt eq ' ' and $tar_pair_nt eq ' ' and 
            $mir_pair_nt eq ' ' and $mir_loop_nt ne ' ') {
            $mispair ++;
            return -1 if $mispair > $mispair_cutoff;
            
            # output alignment symbol
            $tar .= '-';
            $aln_symb .= ' ';
            $mir .= $mir_loop_nt;
            $pos_count ++;
            
            #score
            if ($pos_count >= 2 and $pos_count <= 12) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3_g5 = 1 if $pos_count >= 3 and $pos_count <= 5;
            
            #seed mispairs, cleavage
            if ($pos_count >= 2 and $pos_count <= 7) {
                $seed_mispair ++;
                return -1 if $seed_mispair > $seed_mispair_cutoff;
            } elsif ($pos_count == 8 and defined $mimics) {
                $p8_mispair ++;
            } elsif ($pos_count >= 9 and $pos_count <= 11) {
                $cleavage = 'No';
                $middle_mispair ++ if defined $mimics;
            } elsif ($pos_count >= 12 and $pos_count <= 19 and defined $mimics) {
                $end_mispair ++;
                
                #mimics arm (positions 2-8, 12-19) misp <= $mimics_arm_misp_cutoff
                return -1 if $seed_mispair + $p8_mispair + $end_mispair > $mimics_arm_misp_cutoff;
            }

        } elsif ($tar_loop_nt ne ' ' and $tar_pair_nt eq ' ' and 
                $mir_pair_nt eq ' ' and $mir_loop_nt eq ' ') {
            $mispair ++;
            return -1 if $mispair > $mispair_cutoff;
            $tar .= $tar_loop_nt;
            $aln_symb .= ' ';
            $mir .= '-';
            
            if ($pos_count >= 2 and $pos_count <= 11) { # pos 2-11
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3_g5 = 1 if $pos_count >= 3 and $pos_count <= 4; # pos 2-4
            
            if ($pos_count >= 2 and $pos_count <= 7) {
                $seed_mispair ++;
                return -1 if $seed_mispair > $seed_mispair_cutoff;
            } elsif ($pos_count == 8 and defined $mimics) {
                $p8_mispair ++;
            } elsif ($pos_count >= 9 and $pos_count <= 10) { # pos9-10
                $cleavage = 'No';
                $middle_mispair ++ if defined $mimics;
            } elsif ($pos_count >= 12 and $pos_count <= 19 and defined $mimics) {
                $end_mispair ++;
                return -1 if $seed_mispair + $p8_mispair + $end_mispair > $mimics_arm_misp_cutoff;
            }

        } elsif ($tar_loop_nt ne ' ' and $tar_pair_nt eq ' ' and 
                $mir_pair_nt eq ' ' and $mir_loop_nt ne ' ') {
            $mispair += 1;
            return -1 if $mispair > $mispair_cutoff;
            $tar .= $tar_loop_nt;
            $aln_symb .= ' ';
            $mir .= $mir_loop_nt;
            $pos_count ++;
            
            if ($pos_count >= 2 and $pos_count <= 12) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3_g5 = 1 if $pos_count >= 3 and $pos_count <= 5;
            
            if ($pos_count >= 2 and $pos_count <= 7) {
                $seed_mispair ++;
                return -1 if $seed_mispair > $seed_mispair_cutoff;
            } elsif ($pos_count == 8 and defined $mimics) {
                $p8_mispair ++;
            } elsif ($pos_count >= 9 and $pos_count <= 11) {
                $cleavage = 'No';
                $middle_mispair ++ if defined $mimics;
            } elsif ($pos_count >= 12 and $pos_count <= 19 and defined $mimics) {
                $end_mispair ++;
                return -1 if $seed_mispair + $p8_mispair + $end_mispair > $mimics_arm_misp_cutoff;
            }

        } elsif ( ($tar_pair_nt eq 'U' and $mir_pair_nt eq 'G') or
                ($tar_pair_nt eq 'G' and $mir_pair_nt eq 'U') ) {
            $mispair += 1;
            return -1 if $mispair > $mispair_cutoff;
            $tar .= $tar_pair_nt;
            $aln_symb .= 'o';
            $mir .= $mir_pair_nt;
            $pos_count ++;
            
            if ($pos_count >= 2 and $pos_count <= 12) {
                $score += 1;
            } else {
                $score += 0.5;
            }
            return -1 if $score > $score_cutoff;
            $g3_g5 = 1 if $pos_count >= 3 and $pos_count <= 5;
            
            if ($pos_count >= 2 and $pos_count <= 7) {
                $seed_mispair ++;
                return -1 if $seed_mispair > $seed_mispair_cutoff;
            } elsif ($pos_count == 8 and defined $mimics) {
                $p8_mispair ++;
            } elsif ($pos_count >= 9 and $pos_count <= 11) {
                $cleavage = 'No';
                $middle_mispair ++ if defined $mimics;
            } elsif ($pos_count >= 12 and $pos_count <= 19 and defined $mimics) {
                $end_mispair ++;
                return -1 if $seed_mispair + $p8_mispair + $end_mispair > $mimics_arm_misp_cutoff;
            }

        } else {
            $tar .= $tar_pair_nt;
            $aln_symb .= '|';
            $mir .= $mir_pair_nt;
            $pos_count ++;
        }
    }
    
    #mimics positions 9-11 misp <= 5
    if ( defined $mimics and 
        ($cleavage eq 'Yes' or $middle_mispair > 5 or 
        $seed_mispair + $p8_mispair + $end_mispair > $mimics_arm_misp_cutoff) ) {#reset if needed
        return -1;
    }
        
    my $tarID = $line_arr[0] =~ s/(.+)\_(\d+)$/$1/r;
    my $split_pos = $2;
    my $mirID = $line_arr[2];
    my $rnahyb_pos = $line_arr[6];
    my $rna_start_pos = $split_pos + $rnahyb_pos;

    $tar = reverse $tar;      #change to 5' 3'
    my $mir_r = reverse $mir; #change to 3' 5'
    $aln_symb = reverse $aln_symb;
    
    $score ++ if $g3_g5;
    return -1 if $score > $score_cutoff;
    
    #output detailed info
    open my $rnahyb_details, '>>', "$output_folder/rnahy_out_details" or die $!;
    if ($tab_format) {
        print $rnahyb_details "\#R$rnahy_output_ID\t";
        print $rnahyb_details "$tarID\t$tar\t$mirID\t$mir\t$mispair\t$score\t";
        print $rnahyb_details "$seed_mispair\t$cleavage";
        print $rnahyb_details "\t$rna_start_pos\n";
        
    } else {
        print $rnahyb_details "\#R$rnahy_output_ID", "="x60, "\n";
        printf $rnahyb_details "%-35s\t%s\n", $tarID, "5' $tar 3'";
        printf $rnahyb_details "%-35s\t%s\n", "", "   $aln_symb";
        printf $rnahyb_details "%-35s\t%s\n\n", $mirID, "3' $mir_r 5'";
        printf $rnahyb_details "%-35s\t%s\n", "Total mispairs:", $mispair;
        printf $rnahyb_details "%-35s\t%s\n", "Score:", $score;
        printf $rnahyb_details "%-35s\t%s\n", "Seed mismpairs:",  $seed_mispair;
        printf $rnahyb_details "%-35s\t%s\n", "Cleavage:", $cleavage;
        printf $rnahyb_details "%-35s\t%s\n\n", "Position:", $rna_start_pos;
    }
    close $rnahyb_details;

    return -1 if defined $mimics;

    #pr start-end positions
    my $tar_copy = $tar =~ s/[^\w]//gr;
    my $rna_end_pos = $rna_start_pos + length($tar_copy) - 1;
    my $pr_start_pos = ($rna_start_pos - $rna_start_pos%3) / 3 + 1;
    my $pr_end_pos   = ($rna_end_pos - $rna_end_pos%3) / 3;

    #find algnment positions of pr start-end 
    my ($aln_pr_start, $aln_pr_end, $group);
    for my $k (keys %{$MUSCLE_aln_href}) {
        if ($k =~ /(\S+?)\|$tarID/) {
            ($aln_pr_start, $aln_pr_end) = 
            aln_find_pos($MUSCLE_aln_href->{$k}, $pr_start_pos, $pr_end_pos);
            $group = $1;

            #output rnahy_out file for conservation analysis
            open  my $rnahyb_output, '>>', "$temp_dir/rnahy_out" or die $!;
            print $rnahyb_output "$mirID\t$tarID\t$aln_pr_start\t",
                            "$aln_pr_end\t$group\tR$rnahy_output_ID\n";
            close $rnahyb_output;
        }
    }
}


sub aln_find_pos {
    my ($gene_seq, $start, $end) = @_;

    my $segment1 = substr($gene_seq, 0, $start) =~ s/-//gr;
    my $offset  = $start - length($segment1);
    my $segment2  = substr($gene_seq, $start);
    
    my ($i, $newpos) = (0,0);
    if ($offset > 0) {
        while($segment2 =~ /\w/g) {
            $i++;
            $newpos = pos($segment2);
            last if $i == $offset;
        }
    }
    my $newstart = $start + $newpos;
    
    $offset = $end - $start;
    my $segment3 = substr($gene_seq, $newstart);
    my $j;
    while($segment3 =~ /\w/g) {
        $j++;
        $newpos = pos($segment3);
        last if $j == $offset;
    }
    my $newend = $newstart + $newpos;
    
    return ($newstart, $newend);
}


sub run_FASTA {
    my ($mirName, $mirfile, $cdsfile) = @_;
    
    my $mirAbbr = substr($mirName, 0, 5);
    my %fasta_hash;
    my $find_first = 0;
    my ($geneName, $rna_pair_beg_pos, $rna_pair_end_pos, $rna_start_pos, $rna_end_pos, 
        $mir_end_pos, $mir_rc, $aln, $aln_copy, $tar, $aln_line, $tar_line, $left_pos);
    my (@mir, @tar, @aln);

    STDERR->print ("Running FASTA for $mirName ...  ");
    
    open my $run_fasta, "$fasta_version -A -n -Q -i -U -T $threads -E 100000 $mirfile $cdsfile 1 |" or die $!;
    while (<$run_fasta>) {
        $fasta_output_ID ++;
        chomp;
        unless ($find_first) {
            next unless /^>>(\S+)/;
            $geneName = $1;
            $find_first = 1;
            next;
            
        } else {
            if (/^>>(\S+)/) {
                $geneName = $1;

            } elsif (/Smith-Waterman.+(\d+)-\d+:(\d+)-(\d+)/) {
                $mir_end_pos  = $1;
                $rna_pair_beg_pos = $2; 
                $rna_pair_end_pos = $3;

            } elsif (/^($mirAbbr\-\s+)(\S+)/) {
                $mir_rc = $2;
                
                $left_pos = length($1);
                $aln_line = <$run_fasta>;  #aln symbol
                if ($aln_line =~ /^(\s+)(\S.*\S)(\s*)$/) {
                    $aln = $2;
                    $aln_copy = $aln =~ s/://gr;
                }
                
                $tar_line = <$run_fasta>;
                $tar = substr($tar_line, $left_pos, length($mir_rc));
                $tar =~ s/\s/-/g;    #some targs have unpairings at 5' and 3' end
                
                #primary filtering
                next if (defined $aln_copy and length($aln_copy) > $mispair_cutoff);
                
                FASTA_parser($mirName, $geneName, $mir_rc, $tar, 
                            $rna_pair_beg_pos, $rna_pair_end_pos);
            }
        }
    }
    close $run_fasta;
    STDERR->say ("done.");
}


sub FASTA_parser {
    my ($mirID, $geneID, $mir_rc, $tar, $rna_start_pos, $rna_end_pos) = @_;
    my $mir_rc_copy = $mir_rc =~ s/-//gr;
    return -1 if length($mir_rc_copy) < 10; #remove some low-quality results
    
    (my $mir = reverse uc $mir_rc) =~ tr/AUGC/UACG/;
    $tar = reverse uc $tar;
    my @mir = split //, $mir; #5'-3'
    my @tar = split //, $tar; #3'-5'
    my $aln_symb;
    my $cleavage = 'Yes';

    my ($count, $score, $g3g5, $total_misp, $seed_misp, $p8_misp, $middle_misp, $end_misp)
        = (0,0,0,0,0,0,0,0);

    my %pairing = (
    #class1:pairing, 2:G:U, 3:mismatch, 4:bulge miR strand, 5:bulge tar strand
    'AU' => 'class1',    'UA' => 'class1',    'GC' => 'class1',    'CG' => 'class1',
    'UG' => 'class2',    'GU' => 'class2',    'AG' => 'class3',    'AC' => 'class3',
    'UC' => 'class3',    'GA' => 'class3',    'CA' => 'class3',    'CU' => 'class3',
    'AA' => 'class3',    'UU' => 'class3',    'GG' => 'class3',    'CC' => 'class3', 
    'A-' => 'class4',    'U-' => 'class4',    'G-' => 'class4',    'C-' => 'class4', 
    '-A' => 'class5',    '-U' => 'class5',    '-G' => 'class5',    '-C' => 'class5',
    );
    
    for my $i (0..$#mir) {
        my $mir_nt = $mir[$i];
        my $tar_nt = $tar[$i];
        return -1 unless defined $pairing{"$mir_nt$tar_nt"};
        if ($pairing{"$mir_nt$tar_nt"} eq 'class1') {
            $count ++;
            $aln_symb .= '|';
            
        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class2') {
            $count ++;
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= 'o';
            
            #Score
            if ($count >= 2 and $count <= 12) {
                $score += 1;
             } else {
                $score += 0.5;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 5;
            
            #seed mispairs, cleavage
            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 11) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 12 and $count <= 19 and defined $mimics) {
                $end_misp ++;
                
                #mimics arm (positions 2-8, 12-19) misp <= $mimics_arm_misp_cutoff
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff; 
            }
            
        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class3') {
            $count ++;
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= ' ';
            
            if ($count >= 2 and $count <= 12) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 5;
            
            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 11) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 12 and $count <= 19 and defined $mimics) {
                $end_misp ++;
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff;
            }
            
        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class4') {
            $count ++;
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= ' ';
            
            if ($count >= 2 and $count <= 12) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 5;
            
            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 11) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 12 and $count <= 19 and defined $mimics) {
                $end_misp ++;
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff;
            }

        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class5') {
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= ' ';
            
            if ($count >= 2 and $count <= 11) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 4;
            
            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 10) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 11 and $count <= 19 and defined $mimics) {# pos 12-19
                $end_misp ++;
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff;
            }
        }
    }

    $score ++ if $g3g5;
    return -1 if $score > $score_cutoff;

    #mimics positions 9-11 misp <= 5
    if ( defined $mimics and 
        ($cleavage eq 'Yes' or $middle_misp > 5 or 
        $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff) ) {
        return -1;
    }

    $tar      = reverse $tar; #change to 5'-3'
    my $mir_r = reverse $mir; #change to 3'-5'
    $aln_symb = reverse $aln_symb;

    #output detailed info
    open my $fasta_details, '>>', "$output_folder/fasta_out_details" or die $!;
    if ($tab_format) {
        print $fasta_details "\#F$fasta_output_ID\t";
        print $fasta_details "$geneID\t$tar\t$mirID\t$mir\t$total_misp\t",
                            "$score\t$seed_misp\t$cleavage\t$rna_start_pos\n";
    } else {
        print $fasta_details "\#F$fasta_output_ID", "=="x40, "\n";
        printf $fasta_details "%-35s\t%s\n", $geneID, "5' $tar 3'";
        printf $fasta_details "%-35s\t%s\n", "", "   $aln_symb";
        printf $fasta_details "%-35s\t%s\n\n", $mirID, "3' $mir_r 5'";
        printf $fasta_details "%-35s\t%s\n", "Total mispairs:", $total_misp;
        printf $fasta_details "%-35s\t%s\n", "Score:", $score;
        printf $fasta_details "%-35s\t%s\n", "Seed mispairs:", $seed_misp;
        printf $fasta_details "%-35s\t%s\n", "Cleavage:", $cleavage;
        printf $fasta_details "%-35s\t%s\n\n", "Position:", $rna_start_pos;
    }
    close $fasta_details;
    
    return -1 if defined $mimics;
    
    my $pr_start_pos = ($rna_start_pos - $rna_start_pos%3) / 3 + 1;
    my $pr_end_pos   = ($rna_end_pos - $rna_end_pos%3) / 3;
    
    #output fasta_output file for conservation analysis
    my ($pr_aln_start, $pr_aln_end, $group);
    for my $k (keys %{$MUSCLE_aln_href}) {
        if ($k =~ /(\S+?)\|$geneID$/) { 
            ($pr_aln_start, $pr_aln_end) = 
                aln_find_pos($MUSCLE_aln_href->{$k}, $pr_start_pos, $pr_end_pos);
            $group = $1;
            open  my $fasta_output, '>>', "$temp_dir/fasta_out" or die $!;
            print $fasta_output "$mirID\t$geneID\t$pr_aln_start\t",
                            "$pr_aln_end\t$group\tF$fasta_output_ID\n";
            close $fasta_output;
        }
    }
}


sub conservation {
    my ($infile, $outfile) = @_;
    unless (-e $infile) {
        STDERR->say ("\nNo hits found.\n", "-"x6, " TarHunter analysis completed.");
        exit;
    }
    
    my $infile2 = basename $infile =~ s/(.+)/$1\_details/r;
    system " sort -k 1,1 $infile > $infile\.sorted.tmp ";
    system " join -1 1 -2 2 $infile\.sorted.tmp  $temp_dir/mircluster.sorted.tmp > $infile\.merged.tmp ";
    #merged file format: mirID-geneID-start-end-geneGroup-siteNum-mirCluster 

    my $result_href = parse_merged_file("$infile\.merged.tmp");
    
    #store detailed info into hash
    my %hash_details;
    open my $rnahy_fasta_details, "$output_folder/$infile2" or die $!;
    local $/ = '#';
    my $discarded = <$rnahy_fasta_details>;
    while (<$rnahy_fasta_details>) {
        chomp;
        my ($id, $details);
        if ($tab_format) {
            ($id, $details) = /(.+?)\t(.+)/s or next;
        } else {
            ($id, $details) = /(.+?)\n(.+)/s or next;
        }
        $id =~ s/=//g;
        $hash_details{$id} = $details;
    }
    close $rnahy_fasta_details;
    local $/ = "\n";
    
    #sequenceID as key, speciesID as value
    my %hash_species;
    for my $species (keys %{$species_list_href}) {
        for my $file(glob "$dbs_folder/*.fa") {
            next unless $file =~ /$species\.fa/;
            my $spe_hash = IO::File->new ($file, 'r') or die $!;
            while (<$spe_hash>) {
                chomp;
                $hash_species{$1} = $species if />(\S+)(.*)/;
            }
        }
    }
    
    if ($infile eq "$temp_dir/rnahy_out") {
        STDERR->say ("Outputing RNAhybrid results ...");
    } elsif ($infile eq "$temp_dir/fasta_out") {
        STDERR->say ("Outputing FASTA results ...");
    }
    
    my ($identity_cutoff_pr, $identity_cutoff_nt) = (0.55, 0.65);
    my $final_output = IO::File->new($outfile, 'w') or die $!;
    if ($tab_format) {  #tab format header
        print $final_output "miRID\tmiR sequence\ttargetID\ttarget sequence\tTotal mispair\t",
            "Score\tSeed mispair\tCleavage\tPosition\tSpecies\tIdentity(nt)\%\tIdentity(aa)\%\n";
    }
    
    for my $gene_group (keys %{$result_href}) {
        for my $mir_cluster (keys %{$result_href->{$gene_group}}) { 
            for my $site_num (keys %{$result_href->{$gene_group}{$mir_cluster}}) {
                my $inside_group = $result_href->{$gene_group}{$mir_cluster}{$site_num};
                
                #inside a unique group
                my (@start_arr, @end_arr, %group_species);
                for my $element (keys %{$inside_group}) {
                    my $group_members = $inside_group->{$element};
                    my $species = $hash_species{$group_members->[1]}; #$group_member->[1]:geneID
                    push @start_arr, $group_members->[2];
                    push @end_arr, $group_members->[3];
                    $group_species{$species} ++;
                }
                my $num_of_species = scalar (keys %group_species);
                next if $num_of_species < 2; # at least 2 species
                
                # target site aa boundary
                my $min = min @start_arr;
                my $max = max @end_arr;
                
                my $aa_tar_region;
                open my $site_pr_file, '>>', "$temp_dir/tar_site_pr.fa" or die $!;
                open my $site_nt_file, '>>', "$temp_dir/tar_site_nt.fa" or die $!;
                for my $element (keys %{$inside_group}) {
                    #output target site aa seqs
                    my $gene_id = $inside_group->{$element}->[1];
                    $aa_tar_region = extract_site($min, $max, $gene_id, $gene_group);
                    $aa_tar_region =~ s/-//g;
                    print $site_pr_file ">$gene_id\n$aa_tar_region\n";
                    
                    #output target site nt seqs
                    my (@details_arr, @tar_nt_arr, $tar_nt, $aln_symb);
                    if ($tab_format) {
                        @details_arr = split/\t/, $hash_details{$element};
                        $tar_nt = $details_arr[1]; #site nt seq 5'-3'
                    } else {
                        @details_arr = split/\n/, $hash_details{$element};
                        @tar_nt_arr = split/\s+/, $details_arr[0];
                        $tar_nt = $tar_nt_arr[-2]; #site nt seq 5'-3'
                    }
                    print $site_nt_file ">$gene_id\n$tar_nt\n";
                }
                close $site_pr_file;
                close $site_nt_file;
                
                my ($aln_pr, $site_identity_pr) = aln_identity("$temp_dir/tar_site_pr.fa");
                my ($aln_nt, $site_identity_nt) = aln_identity("$temp_dir/tar_site_nt.fa");
                system "rm $temp_dir/tar_site_nt.fa $temp_dir/tar_site_pr.fa";
                
                next if $site_identity_pr <= $identity_cutoff_pr;
                next if $site_identity_nt <= $identity_cutoff_nt;
                
                my $rounded_identity_nt = sprintf("%.1f", $site_identity_nt * 100);
                my $rounded_identity_pr = sprintf("%.1f", $site_identity_pr * 100);
                
                print $final_output "=="x40, "\n" unless $tab_format;
                
                if ($tab_format) {
                    for my $element (sort {$inside_group->{$a}->[0] cmp $inside_group->{$b}->[0]} 
                                keys %{$inside_group}) {
                        chop $hash_details{$element};
                        print $final_output "$hash_details{$element}\t";
                        print $final_output "$num_of_species\t";
                        print $final_output "$rounded_identity_nt\t$rounded_identity_pr\n";
                    }
                } else {
                    #output details
                    for my $element (sort {$inside_group->{$a}->[0] cmp $inside_group->{$b}->[0]}
                                 keys %{$inside_group}) {
                        print $final_output "$hash_details{$element}\n";
                    }
                    
                    print $final_output "Conserved sites found in $num_of_species species.\n\n";
                    
                    #output target site nt alignment
                    print $final_output "Target site nucleotide acid identity: ", 
                                        $rounded_identity_nt, "\%\n";
                    print $final_output "$aln_nt\n";
                    
                    #output target site aa alignment
                    print $final_output "Target site aa motif identity: ",
                                        $rounded_identity_pr, "\%\n";
                    print $final_output "$aln_pr\n";
                }
            }
        }
    }
}


sub aln_identity {
    my $file = shift;
    my ($identity, $first) = (0,0);
    my ($aln_detail, $len, $symb, %duplicate);
    open my $site_aln_identity, "$muscle_version -in $file -clwstrict -quiet |" or die $!;

    while (<$site_aln_identity>) {
        chomp;
        next if /^CLUSTAL/;
        next if /^\s*$/;
        next if exists $duplicate{$_};
        $duplicate{$_} = 1;
        $aln_detail .= "$_\n"; #target site alignment
        unless ($first) {
            if (/^((\S+?)\s+)\S+/) {
                $len = length $1;
                $first = 1;
            }
        } else {
            next unless /\.|\:|\*/;
            $symb = substr($_, $len);
            my $symb_copy = $symb =~ s/\*//gr;
            $identity = 1-length($symb_copy)/length($symb);
        }
    }
    return ($aln_detail, $identity);
}


sub parse_merged_file {
    my $file = shift;
    my %hash;
    my ($first_in_all, $first_in_group, $siteNo, $aa_shift) = (1, 1, 0, 2);
    my ($curr, $prev);

    open my $merged_file, qq(sort -k 5,5  -k 7,7  -k 3n,3n  -k 4n,4n  $file |); 
    while ($curr = <$merged_file>) {
        $curr =~ s/\r//g;
        chomp $curr;
        next if $curr =~ /^\s+$/;
        if ($first_in_all) {  #first in all
            $prev = $curr;
            $first_in_all = 0;
            next;
        }

        my @a = split/\s+/, $prev;
        my @b = split/\s+/, $curr;
        my ($a2, $a3, $a4, $a5, $a6) = @a[2..6];
        my ($b2, $b3, $b4, $b5, $b6) = @b[2..6];

        if ($a4 eq $b4 and $a6 eq $b6 and abs($a2 - $b2) <= $aa_shift and
            abs($a3 - $b3) <= $aa_shift and $first_in_group ) { #first in group
            $siteNo ++;
            $hash{$a4}{$a6}{$siteNo}{$a5} = [@a[0..3]];  #mirID-geneID-start-end
            $hash{$a4}{$a6}{$siteNo}{$b5} = [@b[0..3]];  
            $first_in_group = 0;

        } elsif ($a4 eq $b4 and $a6 eq $b6 and abs($a2 - $b2) <= $aa_shift and
            abs($a3 - $b3) <= $aa_shift and $first_in_group == 0) {
            $hash{$a4}{$a6}{$siteNo}{$b5} = [@b[0..3]];

        } else {
            $first_in_group = 1;
        }
        $prev = $curr;
    }
    close $merged_file;
    return \%hash;
}


sub extract_site {
    my ($min, $max, $gene, $group) = @_;
    my ($tarsite_region, $region, $reg1, $reg2);
    my $seq = $MUSCLE_aln_href->{"$group\|$gene"};

    # target site aa motif
    $tarsite_region = substr($seq, $min - 1, $max - $min + 1);
    return $tarsite_region;
}


sub remove_duplicate {
    my ($infile, $outfile) = @_;
    my %result;
    my $delimiter = "=="x40;
    local $/ = "$delimiter";
    my $infile_fh = IO::File->new ($infile, 'r') or die $!; 
    while (<$infile_fh>) {
        chomp;
        next if /^\s*$/;
        $result{$_} = 1;
    }
    $infile_fh->close;
    
    my $outfile_fh = IO::File->new ($outfile, 'w') or die $!;
    for my $group (keys %result) {
        print $outfile_fh "$delimiter\n$group";
    }
}


sub mimics_nctarg_conservation {

    local $/ = '#';
    my %mir_group;
    my $mir_cluster_file;
    if ($nonmirbase_file) { $mir_cluster_file = IO::File->new ("$temp_dir/cdhit/mir_cluster.txt") }
    else                  { $mir_cluster_file = IO::File->new ("$dir/miRs/mir_cluster.txt", 'r') }
    my $omitted = <$mir_cluster_file>;
    while(<$mir_cluster_file>) {
        chomp;
        s/\r//g;
        my @arr = split/\n/, $_;
        for my $member (@arr[1..$#arr]) {
            next if $member =~ /^\s*$/;
            $member =~ s/^\s*|\s*$//;
            $mir_group{$member} = $arr[0];
        }
    }
    $mir_cluster_file->close;

    # store fasta/rnahy output details into hash
    my ($delimiter, $records, $fasta_rnahy_file, $final_out);
    my (%mimics_species, %info, %hash_details);
    
    if ($program_option == 1) {
        $fasta_rnahy_file = "$output_folder/fasta_out_details";
        $final_out = "$output_folder/fastaout_nc_conserv";
    } elsif ($program_option == 2) {
        $fasta_rnahy_file = "$output_folder/rnahy_out_details";
        $final_out = "$output_folder/rnahyout_nc_details";
    }
    
    unless (-e $fasta_rnahy_file) {
        STDERR->say ("\nNo hits found.\n", "-"x6, " TarHunter analysis completed.");
        exit;
    }
    
    $tab_format ? ($delimiter = "\t") : ($delimiter = "\n");
    my $details_outfile = IO::File->new($fasta_rnahy_file, 'r') or die $!;
    while (<$details_outfile>) {
        s/\r//g;
        chomp;
        next if /^\s*$/;
        my @arr = split/$delimiter/, $_;
        my $ID = $arr[0] =~ s/=//gr;
        my $geneID = $1 if $arr[1] =~ /^(\S+).*/;
        my $mirID = $1 if $arr[3] =~ /^(\S+).*/;
        my $start_pos = $1 if $arr[9] =~ /.*?(\d+)$/;
    
        $mimics_species{$1} = 1 if $mirID =~ /^(\w+)-miR/;
    
        $info{$ID} = [$mirID, $geneID, $start_pos];
    
        my @arr2 = (split/$delimiter/, $_, 2) or next;
        $hash_details{$ID} = $arr2[1];
    }
    $details_outfile->close;
    local $/ = "\n";

    # get all sequences from mimics result
    my %mimics_all_seq;
    for my $species (keys %mimics_species) {
        my $seq_href = read_seq_file("$dbs_folder/$species.fa");
        @mimics_all_seq{keys %{$seq_href}} = values %{$seq_href};
    }

    # extract upstream and downstream 50-bp (default) seqs
    for my $ID (sort keys %info) {
        my $mirID = $info{$ID}->[0];
        my $geneID = $info{$ID}->[1];
        my $pos = $info{$ID}->[2];
 
        my $seq = $mimics_all_seq{$geneID};
        my $upstream = ($pos > int($nc_len/2)) ? ($pos - int($nc_len/2)) : 1;
        my $newseq = substr($seq, $upstream - 1, $nc_len);

        my $species = $1 if $mirID =~ /^(\w+)-miR/;

        open my $nc_seq, '>>', "$temp_dir/mimics/$species.fa";
        print $nc_seq ">$ID\|$mirID\|$geneID\n$newseq\n";
        close $nc_seq;
    }

    # run usearch
    my %checked;
    my $db_list = IO::File->new ("$temp_dir/mimics/db_list.txt", 'w') or die $!;
    for my $species1 (keys %mimics_species) {
        my $db1 = "$temp_dir/mimics/$species1.fa";
        next unless (-e $db1);
        system "$usearch_version -makeudb_usearch $db1 -output $temp_dir/mimics/$species1.udb >>$temp_dir/msg.txt 2>&1";
            
        for my $species2 (keys %mimics_species) {
            next if $species1 eq $species2;
            next if exists $checked{"$species2$species1"};
            my $db2 = "$temp_dir/mimics/$species2.fa";
            next unless (-e $db2);
            $checked{"$species1$species2"} = 1;
            
            $db_list->print ("$db2\t$temp_dir/mimics/$species1.udb\t$temp_dir/mimics/$species2\_$species1\_result\n");
        }
    }
    $db_list->close;
    
    system "parallel --no-notice --jobs $jobs --colsep '\t' -a $temp_dir/mimics/db_list.txt $usearch_version -usearch_global {1} -db {2} -strand plus -threads $threads -id $iden_cutoff -blast6out {3} >>$temp_dir/msg.txt 2>&1";

    my @result_files = glob "$temp_dir/mimics/*result";
    STDERR->say ("\n", "-"x6, " TarHunter analysis completed.") and exit if @result_files < 1;
    system "cut -f1,2 $temp_dir/mimics/*result | sort -u > $temp_dir/mimics/ublast_pairs";


    # parsing usearch result
    my %pair;
    my $ublast_pair_file = IO::File->new("$temp_dir/mimics/ublast_pairs", 'r');
    while (<$ublast_pair_file>) {
        chomp;
        my @line = split;
        my @id1 = split/\|/, $line[0];
        my @id2 = split/\|/, $line[1];

        if ( $mir_group{$id1[1]} eq $mir_group{$id2[1]} ) {
            $pair{$id1[0]}{$id2[0]} = 1;
            $pair{$id2[0]}{$id1[0]} = 1;
        }
    }
    $ublast_pair_file->close;

    my ($elements_href, $group_num, $group_href);
    while ( my ($memb, $v) = each %pair ) {
        my @elems = keys %{ recursive_grouping(\%pair, $memb, $elements_href) };
        $group_num ++;

        delete @pair{@elems};
        for (@elems) {
            $group_href->{$group_num}{$_} = 1;
        }
    }

    my $output = IO::File->new ($final_out, 'w');
    for my $groupID (keys %{$group_href}) {
        print $output ("=="x40, "\n") unless $tab_format;
        my $seqs = IO::File->new ("$temp_dir/mimics/seqs.fa", 'w');
        
        for my $ID ( keys %{$group_href->{$groupID}} ) {
            print $output "$hash_details{$ID}\n";
            unless ($tab_format) {
                my $gene = $info{$ID}->[1];
                my $pos = $info{$ID}->[2];
                my $seq = $mimics_all_seq{$gene};
                my $upstream = ( $pos > int($nc_len/2) ) ? ($pos - int($nc_len/2)) : 1;
                my $seq_100bp = substr($seq, $upstream - 1, $nc_len);
                print $seqs ">$gene\n$seq_100bp\n"; 
            }
        }
        
        unless ($tab_format) {
            open my $seqs_aln, "$muscle_version -in $temp_dir/mimics/seqs.fa -clwstrict -quiet |" or die $!;
            print $output "Alignment of target sites:\n";
            while (<$seqs_aln>) {
                chomp;
                next if /^CLUSTAL/;
                print $output "$_\n";
            }
        }
    }
}


sub recursive_grouping {
    my ($gene_pair, $gene, $elements_href) = @_;
    
    if (exists $elements_href->{$gene}) {
        return;
    } else {
        $elements_href->{$gene} = 1;
    }
    
    for my $k (keys %{$gene_pair->{$gene}}) {
        recursive_grouping($gene_pair, $k, $elements_href);
    }
    return $elements_href;
}

__END__