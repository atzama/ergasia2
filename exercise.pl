#!/usr/bin/perl
use strict;
use warnings;

my $dna = "ATGAAATGAGGGTAGTAAATGTGA";

# Forward strand
print ">>> Forward strand:\n";
print "Sequence: $dna\n";
find_orfs($dna);

# Reverse complement strand
my $revcomp = reverse_complement($dna);
print "\n>>> Reverse complement strand:\n";
print "Sequence: $revcomp\n";
find_orfs($revcomp);

# Function to find ORFs and translate them into proteins
sub find_orfs {
    my ($seq) = @_;
    while ($seq =~ /(ATG(?:...)*?(?:TAA|TAG|TGA))/g) {
        my $orf = $1;
        print "ORF: $orf\n";
        my $protein = translate($orf);
        print "Protein: $protein\n";
    }
}   

# Function that calculates the reverse complement of a DNA sequence
sub reverse_complement {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;    # Complement
    return reverse($seq);             # Reverse
}

# Function that translates a DNA ORF into a protein sequence
sub translate {
    my ($orf) = @_;
    my %table = (
        'ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M',
        'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
        'AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K',
        'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R',
        'CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
        'CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
        'CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q',
        'CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R',
        'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V',
        'GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
        'GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E',
        'GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G',
        'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S',
        'TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L',
        'TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGA'=>'*'
    );
    my $protein = '';
    for (my $i = 0; $i < length($orf) - 2; $i += 3) {
        my $codon = substr($orf, $i, 3);
        $protein .= $table{$codon} || 'X';  # Use 'X' if codon is unknown
    }
    return $protein;
}
