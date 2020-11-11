# Copyright 2020 EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# This script takes a set of species input dbs on a server, cuts their genomes up into fixed length regions and then encodes the exon/intron/cds features for each protein coding gene in each sampled region. As output you get one file with each sampled region and another with the encoded version of them.
#
# At the moment it only looks at genes on the forward strand, as there are a bunch of extra bits of code and testing needed to correctly encode the reverse strand. But we can look at adding that in later (which will effectively double our possible training set).
#
# The way I have been thinking about it is on a per base labelling, i.e. can we use the context of the entire sequence to decide what an individual base should be labelled. So the encoded sequence is really just the set of labels for each base, but the labels can’t be worked out without the surrounding genomic sequence (hence the large genomic regions).
#
# I have given a list of potential high quality dbs to train on in the script itself. I did a run last night on them and the output is here:
# /hps/nobackup2/production/ensembl/fergal/coding/gene_finder/run_1_genomic_seqs.fa
# /hps/nobackup2/production/ensembl/fergal/coding/gene_finder/run_1_encoded_seqs.fa
#
# As the sequences are 500kb long each it’s a little hard to visualise. For the encoded sequences there’s naturally large sections of Os (representing intergenic region), with the gene features making up small regions. Actually in the file listed above you can see a bunch of single exon genes encoded in arabidopsis at the very start of the file “OOOOSSSCCCCCCCCC….CCCCCCCCEEEOOOOO”. Most genes in general have introns, but you will still see plenty of single exon ones too.
#
# There’s also some code at the bottom of the script that you can enable if you want to print out the raw and encoded exons, which makes things easier to visualise.


use warnings;
use strict;
use feature 'say';
use List::Util 'shuffle';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);


my $coord_system = 'toplevel';

my $user = 'ensro';
my $host = 'mysql-ens-vertannot-staging';
my $port = 4573;
my $pass;
my $db_input_file; # The input file containing all the database names to get sampled on the host server
my $slice_target_length = 500000; # The target length of the region to be analysed for genes
my $slice_overlap = 250000; # The amount of overlap between regions (to capture genes cut by regions boundaries)
my $random_sample = 1; # If you want to randomly sub sample regions from each genome
my $num_samples = 2000; # The max number of random region samples per genome, if randomly sub sampling
my $genomic_output_file = "genomic_seqs.fa"; # Where the all the genomic sequences for each sampled region are written
my $encoded_output_file = "encoded_seqs.fa"; # Where the encoded version of the genomic sequences are written

my $options = GetOptions ("user|dbuser|u=s"       => \$user,
                          "host|dbhost|h=s"       => \$host,
                          "port|dbport|P=i"       => \$port,
                          "db_input_file=s"       => \$db_input_file,
                          "slice_length=i"        => \$slice_target_length,
                          "slice_overlap=i"       => \$slice_overlap,
                          "random_sample!"        => \$random_sample,
                          "num_samples=i"         => \$num_samples,
                          "genomic_output_file=s" => \$genomic_output_file,
                          "encoded_output_file=s" => \$encoded_output_file);


###############################
# Example input file contents
###############################
#arabidopsis_thaliana_core_50_103_11
#caenorhabditis_elegans_core_103_269
#drosophila_melanogaster_core_103_9
#homo_sapiens_core_103_38
#mus_musculus_core_103_39
#saccharomyces_cerevisiae_core_103_4
#triticum_aestivum_core_50_103_4
#zea_mays_core_50_103_7

unless($db_input_file && -e $db_input_file) {
  die "Input file does not exist";
}

open(IN,$db_input_file);
my @dbs = <IN>;
close IN;

open(GENOMIC_OUT,">".$genomic_output_file);
open(ENCODED_OUT,">".$encoded_output_file);


#my @dbs = ('mus_musculus_core_103_39');

my $seq_counter = 1; # Counter for sequence headers in output file, equivalent genomic and encoded sequences will have the same header

# Loop through the input dbs
foreach my $dbname (@dbs) {
  chomp($dbname);
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -port    => $port,
    -user    => $user,
    -host    => $host,
    -dbname  => $dbname);

  # Get all slices
  my $slice_adaptor = $db->get_SliceAdaptor();
  my $slices = $slice_adaptor->fetch_all('toplevel');

  my $gene_adaptor = $db->get_GeneAdaptor();

  my $sub_slices = [];
  foreach my $slice (@$slices) {
    # Just for testing, this only looks at chromosome 4 (assuming your test has a chr 4)
    # unless($slice->seq_region_name =~ /^\d+$/ && $slice->seq_region_name eq '4') {
    #   next;
    # }

    # Cut the toplevel slices up into sub slices, based on a target, length and overlap
    push(@$sub_slices,@{create_sub_slices($slice,$slice_adaptor,$slice_target_length,$slice_overlap,1)});

    # Some code for creating reverse versions of the forward strands
    # Could be a basis for encoding the genes on the reverse strand, however more code is needed to ensure the feature coordinates
    # are also correctly reversed
    # my $reverse_slice_name = $slice->name();
    # chop($reverse_slice_name);
    # $reverse_slice_name .= '-1';
    # my $reverse_slice = $slice_adaptor->fetch_by_name($reverse_slice_name);
    # push(@$sub_slices,@{create_sub_slices($reverse_slice,$slice_adaptor,$slice_target_length,$slice_overlap,-1)})
  } # End foreach my $slice


  # If you elect to randomly sub sample, this will shuffle the indices and then take the equivalent number of shuffled slices
  # If there are less sub slices than the requested number of samples, it just takes all sub slices
  if($random_sample && $num_samples < scalar(@$sub_slices)) {
    my $random_sub_slices = [];
    my @shuffled_indexes = shuffle(0..$#$sub_slices);
    my @pick_indexes = @shuffled_indexes[ 0 .. $num_samples - 1 ];

    foreach my $index (@pick_indexes) {
      push(@$random_sub_slices,${$sub_slices}[$index]);
    }

    $sub_slices = $random_sub_slices;
  }

  # Once you have whatever set of sub slices for the current genome, process and encode the genomic sequence
  # Write the genomic and encoded sequences to the output files
  foreach my $sub_slice (@$sub_slices) {
    my ($genomic_seq,$encoded_seq) = process_slice($sub_slice,$slice_target_length,$gene_adaptor);
    say GENOMIC_OUT ">".$seq_counter."\n".$genomic_seq;
    say ENCODED_OUT ">".$seq_counter."\n".$encoded_seq;
    $seq_counter++;
  }

} # End foreach my $dbname

close GENOMIC_OUT;
close ENCODED_OUT;

exit;


sub create_sub_slices {
  my ($slice,$slice_adaptor,$slice_target_length,$slice_overlap,$strand) = @_;

  # Takes a target slice length and overlap and then cuts the original slice up based on this
  # Later on sub slices that are shorter than the target length are padded up to the target length

  # If the slice is already smaller than the target length, then just return it
  if($slice->length < $slice_target_length) {
    return([$slice]);
  }

  # Cut up the slice into sub slices
  my $sub_slices = [];
  my $start = 1;
  my $end = $slice_target_length;
  while($end <= $slice->length()) {
    my $sub_slice = $slice_adaptor->fetch_by_region( 'toplevel',$slice->seq_region_name,$start,$end,$strand);
    push(@$sub_slices,$sub_slice);
    $start = $end - $slice_overlap + 1;
    $end += $slice_target_length - $slice_overlap;
  }

  if($end > $slice->length()) {
    $end = $slice->length();
    my $sub_slice = $slice_adaptor->fetch_by_region( 'toplevel',$slice->seq_region_name,$start,$end,$strand);
    push(@$sub_slices,$sub_slice);
  }

  return($sub_slices);
}



sub process_slice {
  my ($slice,$slice_target_length,$gene_adaptor) = @_;

  # Something like this can be used for testing a particular slice
  # unless($slice->name eq 'chromosome:GRCm39:4:4500001:5000000:1') {
  #   return;
  # }

  # This takes the current slice and processes it. If it's less than the target length, the underlying sequence gets padded
  # to ensure all regions are equivalent length
  my $genomic_seq = $slice->seq();
  my $length_diff = $slice_target_length - $slice->length();
  my $padding = "T" x $length_diff;
  if($padding) {
    $genomic_seq .= $padding;
  }

  # Get all the genes from the slice
  my $genes = $slice->get_all_Genes();
  say "Processing: ".$slice->name." (".scalar(@$genes)." genes)";

  my $canonical_transcripts = [];
  foreach my $gene (@$genes) {
    my $biotype = $gene->biotype;

    # Note, at the moment it's skipping genes on the reverse strand. These should be encoded on the reverse complement sequence
    # So instead of just having one slice per region, have two slices per region, one forward and one reverse
    # As genes can overlap on opposite strands (and thus a base could in in two genes, one on each strand), it's better to
    # process the strands independently. This will also help with learning, instead of having to learn that genes can be in either
    # orientation, we will always interpret in a 5' to 3' direction
    unless($biotype eq 'protein_coding' && $gene->strand == 1) {
      next;
    }

    # At this point we know the current gene is coding and on the forward strand
    # First we check if the gene goes off the edge of the slice, if so we mask the gene out on the region and move onto the next gene
    $genomic_seq = mask_seq($gene,$genomic_seq,$slice);
    if($gene->seq_region_start < $slice->seq_region_start() || $gene->seq_region_end() > $slice->seq_region_end()) {
      next;
    }

    # At this point we know the gene is completely within the region, so we get the canonical transcript. The canonical transcript
    # becomes the basis for encoding the features for that gene on the encode version of the genomic region
    # For the moment, we will just collect the canonical transcripts for all the genes that get to this point into an array
    # and then process them outside the loop
    my $canonical_transcript = $gene->canonical_transcript();
    unless($canonical_transcript) {
      die "Didn't find a canonical";
    }

    push(@$canonical_transcripts,$canonical_transcript);
  } # End foreach my $gene

  # We now have all the canonical transcripts for genes that are fully within the region. We use their features to now create the
  # encoded sequence
  my $encoded_seq = encode_seq($genomic_seq,$canonical_transcripts,$slice);

  # Return back both the genomic sequence and the encoded version of it. Note that the genomic sequence may also have been altered
  # in that it may have been padded to bring it up to the target length and also edge regions where the genes went off the edge of
  # the region will have been masked out
  return($genomic_seq,$encoded_seq);
}


sub mask_seq {
  my ($gene,$genomic_seq,$slice) = @_;

  # If the gene is over the edge of the sub slice, mask the sequence out so that it doesn't confuse the training
  # The problem with leaving them in would be that we would have part of the signature of the gene in the region
  # but not all of it, so even if we didn't encode it, it would still have signal present in the genomic sequence
  my $seq_region_start = $slice->seq_region_start();
  my $start = $gene->seq_region_start() - $seq_region_start;
  my $end = $gene->seq_region_end() - $seq_region_start;

  if($start < 0 && $end > $slice->length()) {
    return("T" x length($genomic_seq));
  }

  if($start < 0) {
    substr($genomic_seq,0,$end) = "T" x $end;
  }

  if($end > $slice->length()) {
    substr($genomic_seq,$start,length($genomic_seq) - $start) = "T" x (length($genomic_seq) - $start);
  }

  return($genomic_seq);
}


sub encode_seq {
  my ($genomic_seq,$canonical_transcripts,$slice) = @_;

  # This takes the genomic sequence and encodes it based on the following encoding:
  # O - Bases not representing part of a gene (specifically not representing part of the CDS and introns)
  # C - Bases in the coding region of a gene
  # S - Bases part of the start codon of the coding region (three such consective bases per gene)
  # E - Bases part of the stop codon of the coding region (three such consective bases per gene)
  # I - Bases part of the introns of the gene
  # J - Bases part of the splice junction of an intron, there should be a junction on each end of the intron

  # These are genomic coords, as the API is a little inconsistent sometimes in terms of giving slice coords
  # versus genomic coords
  my $seq_region_start = $slice->seq_region_start();
  my $seq_region_end = $slice->seq_region_end();


  # Start by encoding the sequence as Os
  my $encoded_seq = "O" x length($genomic_seq);

  # Then loop through the transcripts and encode the features on top
  foreach my $transcript (@$canonical_transcripts) {
    my $translateable_seq = $transcript->translateable_seq;
    unless($translateable_seq =~ /^ATG/ && ($translateable_seq =~ /TAA$/ ||
                                            $translateable_seq =~ /TGA$/ ||
                                            $translateable_seq =~ /TAG$/)) {
      # It the gene looks unusual (i.e, no start/end) effectively mask it out in both the genomic and encoded seq
      my $start = $transcript->start() - $seq_region_start;
      substr($genomic_seq,$start,$transcript->length()) = "T" x $transcript->length();
      substr($encoded_seq,$start,$transcript->length()) = "O" x $transcript->length();
      next;
    }

    my @cds_exons = @{$transcript->get_all_CDS()};

    for my $i (0..$#cds_exons) {
      my $exon = $cds_exons[$i];
      my $start = $exon->start() - $seq_region_start;
      my $end = $exon->end() - $seq_region_end;

      my $coded_exon = "C" x $exon->length();
      substr($encoded_seq,$start,length($coded_exon)) = $coded_exon;
    }

    my $coding_start = $transcript->coding_region_start() - $seq_region_start;
    my $coding_end = $transcript->coding_region_end() - $seq_region_start;
    substr($encoded_seq,$coding_start,3) = "SSS";
    substr($encoded_seq,$coding_end-2,3) = "EEE";

    my $introns = $transcript->get_all_Introns();
    foreach my $intron (@$introns) {
      unless($intron->next_Exon->is_coding($transcript) && $intron->prev_Exon->is_coding($transcript)) {
        next;
      }

      my $start = $intron->start() - $seq_region_start;
      my $end = $intron->end() - $seq_region_start;
      my $coded_intron = "I" x $intron->length();
      substr($encoded_seq,$start,length($coded_intron)) = $coded_intron;
      substr($encoded_seq,$start,3) = "JJJ";
      substr($encoded_seq,$end-2,3) = "JJJ";
    }
  }

  # Check encoding, if you enable this code it will show the exon genomic seqs and the exon encoded seq
  # Checking this is a relatively easy to see the coding region and the start/stop codons are being encoded
  # correctly. Easier than looking at the full seq as introns are quite long
  # foreach my $transcript (@$canonical_transcripts) {
  #   foreach my $exon (@{$transcript->get_all_Exons()}) {
  #     unless($exon->is_coding($transcript)) {
  #       next;
  #     }
  #
  #     my $start = $exon->start() - $seq_region_start;
  #     my $exon_raw = substr($genomic_seq,$start,$exon->length());
  #     my $exon_encoded = substr($encoded_seq,$start,$exon->length());
  #     say "R: ".$exon_raw;
  #     say "E: ".$exon_encoded;
  #   }
  # }

  return($encoded_seq);
}
