#!/usr/bin/env perl

# generate_consensus_from_fasta_and_vcf.pl
#
# Takes reference sequences in FASTA format and a VCF file containing SNPs
# (but not indels) and replaces the reference bases with alternate bases
# For Heliconius speciation
# Author: John Davey john.davey@ed.ac.uk
# Begun 7/4/11

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;

# Autoflush output so reporting on progress works
$| = 1;

my $ref_filename = "";
my $vcf_filename = "";
my $indels;
my $options_okay = GetOptions(
    'ref=s'  => \$ref_filename,
    'vcf=s'  => \$vcf_filename,
    'indels' => \$indels,
);

croak
"\nUsage: perl generate_consensus_from_fasta_and_vcf.pl -r reference -v vcf -i process_indels\n"
  if !$options_okay;

open my $ref_file, '<', $ref_filename
  or croak "Can't open $ref_filename: $OS_ERROR!\n";

open my $vcf_file, '<', $vcf_filename
  or croak "Can't open $vcf_filename: $OS_ERROR!\n";

my $first_vcf_line = "";
while ( my $vcf_line = <$vcf_file> ) {
    if ( $vcf_line !~ /^#/ ) {
        $first_vcf_line = $vcf_line;
        last;
    }
}

my $scaffold;
my @scf_lines;
while ( my $ref_line = <$ref_file> ) {

    if ( $ref_line =~ /^>(.+)/ ) {    # header line

        # If not first header line, output the scaffold
        if ( @scf_lines != 0 ) {
            my $new_scaffold = $1;

            # Get variants, if there are any (if the already loaded
            # VCF line is for this scaffold)
            my %variants;
            if ( $first_vcf_line =~ /^$scaffold/ ) {
                $first_vcf_line =
                  load_vcf_lines( \%variants, $vcf_file, $scaffold,
                    $first_vcf_line, $indels );
            }

            # Output scaffold incorporating any variants found
            process_scaffold( \@scf_lines, \%variants, $scaffold );

            # Initialise new scaffold and output new scaffold header
            $scaffold  = $new_scaffold;
            @scf_lines = ();
            print $ref_line;
        }
        else {    # First scaffold, just print header
            $scaffold = $1;
            print $ref_line;
        }
    }
    else {        # Sequence line; add to sequence lines for this scaffold
        chomp $ref_line;
        push @scf_lines, $ref_line;
    }
}

# Process last scaffold
my %variants;
if ( $first_vcf_line = /^$scaffold/ ) {
    $first_vcf_line = load_vcf_lines( \%variants, $vcf_filename, $scaffold, $first_vcf_line,
        $indels );
}

process_scaffold( \@scf_lines, \%variants, $scaffold );
if ($first_vcf_line ne "") {print STDERR "VCF file not empty!\n";}

close $ref_file;
close $vcf_file;
exit;

sub load_vcf_lines {
    my ( $variants_ref, $vcf_file, $scaffold, $first_vcf_line, $indels ) = @_;

    process_vcf_line( $variants_ref, $first_vcf_line, $indels );

    while ( my $vcf_line = <$vcf_file> ) {
        if ( $vcf_line =~ /^$scaffold/ ) {
            process_vcf_line( $variants_ref, $vcf_line, $indels );
        }
        else {
            return $vcf_line;
        }
    }

    return "";
}

sub process_vcf_line {
    my ( $variants_ref, $vcf_line, $indels ) = @_;

    chomp $vcf_line;
    my @vcf_fields = split /\t/, $vcf_line;
    my $alternate  = $vcf_fields[4];
    my $attributes = $vcf_fields[7];

    # If not loading indels, skip them
    return if ( ( !$indels ) && ( $attributes =~ /^INDEL/ ) );

    # If no alternate and not an indel (ie reference call), skip
    return
      if ( ( $attributes !~ /^INDEL/ )
        && (( $alternate eq "." ) || ($alternate eq "N")) );

    # Only take the first alternate
    $alternate = ( split /,/, $alternate )[0];

    my $position  = $vcf_fields[1];
    my $reference = $vcf_fields[3];

    $variants_ref->{$position}{ref} = $reference;
    $variants_ref->{$position}{alt} = $alternate;

    return;
}

sub process_scaffold {
    my ( $scf_lines_ref, $variants_ref, $scaffold ) = @_;

    my $pos = 0;

    # Process and output scaffold
    foreach my $scf_line ( @{$scf_lines_ref} ) {

        my @bases = split //, $scf_line;

        for my $base (@bases) {
            $pos++;
            if ( defined $variants_ref->{$pos} ) {
                my $var = $variants_ref->{$pos};

                # Preserve ref to current pos
                if (   ( length( $var->{ref} ) == 1 )
                    && ( $var->{ref} ne $base ) )
                {
                    print STDERR "$scaffold:$pos $base is not $var->{ref}\n";
                }

                if ( $var->{alt} ne "." ) {
                    print $var->{alt};
                }

                # Remove indel bases, accepting indels called within
                # the range of the current indel
                my $skipbases = length( $var->{ref} ) - 1;
                while ( $skipbases > 0 ) {
                    shift @bases;
                    $pos++;
                    $skipbases--;

                    # If the new position is itself an indel,
                    # and the new indel is longer than the existing indel,
                    # add its length to the number of bases to skip
                    if ( defined $variants_ref->{$scaffold}{$pos} ) {
                        if (
                            length( $variants_ref->{$scaffold}{$pos}{ref} ) >
                            length( $var->{ref} ) )    # New pos v old pos
                        {
                            $skipbases +=
                              length( $variants_ref->{$scaffold}{$pos}{ref} ) - 1;
                        }
                    }
                }
            }
            else {
                print $base;
            }
        }
        print "\n";
    }

}
