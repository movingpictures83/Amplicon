use strict;
use Bio::SeqIO;

my $s;
my $in;
my $st;
my $fout;
my $lout;
my $forPrimer;
my $revPrimer;
sub cleanPrimer {
    # fix primers so that they have no IUPAC codes
    my ($p) = @_;
    # print "cleanPrimer: $p:\n";
    # convert to upper case
    $p =~ tr/acgtuUrymkwsbdhvn/ACGTTTRYMKWSBDHVN/;
    $p =~ s/R/\[AG\]/g;
    $p =~ s/Y/\[CT\]/g;
    $p =~ s/M/\[CA\]/g;
    $p =~ s/K/\[TG\]/g;
    $p =~ s/W/\[TA\]/g;
    $p =~ s/S/\[CG\]/g;
    $p =~ s/B/\[CTG\]/g;
    $p =~ s/D/\[ATG\]/g;
    $p =~ s/H/\[ATC\]/g;
    $p =~ s/V/\[ACG\]/g;
    $p =~ s/N/\[ACGT\]/g;
    # print "cleanPrimer: $p;\n";
    return $p;
}


sub input {
   my $inputfile;
   my $fin;
   open(my $inputfile, '<:encoding(UTF-8)', @_[0])
     or die "Could not open file '@_[0]' $!";
   $fin = <$inputfile>;
   chomp($fin);
   $fout = <$inputfile>;
   chomp($fout);
   $lout = <$inputfile>;
   chomp($lout);
   $forPrimer = <$inputfile>;
   chomp($forPrimer);
   $revPrimer = <$inputfile>;
   chomp($revPrimer);
   $in = Bio::SeqIO->newFh(-file => "<$fin", -format => 'fasta');
   return;
}

sub run {
   return;
}

sub output {
open(FH, '>', @_[0]) or die $!;
my $out = Bio::SeqIO->new(-file => ">$fout", -format => 'fasta');
unless (open(LFILE, ">$lout")) {
    print "File $lout does not exist"; exit;
}

$forPrimer = cleanPrimer($forPrimer);
$revPrimer = cleanPrimer($revPrimer);
# find reverse complement of <RevPrimer> 
my $revPrimer2 = reverse($revPrimer);
$revPrimer2 =~ tr/ACGT\[\]/TGCA\]\[/;

# find lengths of primers
my $lF = length($forPrimer);
my $rF = length($revPrimer);

# some output messages sent to screen
print FH "Primers: \nFORWARD:\t$forPrimer ($lF)\nREVERSE:\t$revPrimer ($rF)\n";
print FH "REVERSE:\t$revPrimer2 (complemented)\n";
print FH "FOUND IN FOLLOWING SEQUENCES:\n";
my $seq;
my $count = 0;

print FH "Length\tDomain\tGenus\tSpecies\tStrain\tGene ID\tNumber\tFindCode\n";
while ($seq = <$in>) {
    $count++;
    my $d = $seq->desc();
    my $s = $seq->seq();

    # First processing the description of the fasta sequence
    my $domain = "";
    if ($d =~ "archae") {
	$domain = "Arch ";
    } elsif ($d =~ "prokaryote") {
	$domain = "? ";
    } else {
	$domain = "B ";
    }
    my $genus = "-";
    my $species = "-";
    my $strain = "-";
    my $gid = "-";

    # Just a kludgy fix
    if ($d =~ "unidentified;") {
	$d =~ s/unidentified;/unidentified bacterium;/;
    }

    my @DA = split(/;/, $d);
    if ($#DA >= 2) {
	$strain = $DA[1];
	$strain =~ s/ //g; # Do I really need this?
	$gid = $DA[$#DA];
	my @parts = split(/ /, $DA[0]);
	if ($#parts >= 0) {
	    if (($parts[0] ne "uncultured") && ($parts[0] ne "unidentified")) {
		$genus = $parts[0];
		$species = $parts[1];
	    } elsif (($parts[1] ne "bacterium") && ($parts[1] ne "archaeon") &&
		     ($parts[1] ne "eubacterium") && ($parts[1] ne "prokaryote")) {
		$genus = $parts[1];
	    }
	}
    } elsif ($#DA == 1) {
	$gid = $DA[1];
	my @parts = split(/ /, $DA[0]);
	if ($#parts >= 0) {
	    if (($parts[0] ne "uncultured") && ($parts[0] ne "unidentified")) {
		$genus = $parts[0];
		if ($#parts >= 1) {
		    if (($parts[1] ne "sp.") && ($parts[1] ne "str.")
			&& ($parts[1] ne "bacterium")) {
			$species = $parts[1];
			if ($#parts >= 2) {
			    $strain = $parts[2];
			    if ($#parts >= 3) {
				$strain = $strain.$parts[3]; # what if there is more
			    }
			}
		    } else {
			if ($#parts >= 2) {
			    $strain = $parts[2];
			    if ($#parts >= 3) {
				$strain = $strain.$parts[3]; # what if there is more
			    }
			}
		    }			
		}
	    } elsif (($parts[1] ne "bacterium") && ($parts[1] ne "archaeon") &&
		     ($parts[1] ne "eubacterium") && ($parts[1] ne "prokaryote")) {
		$genus = $parts[1];
		if ($#parts >= 2) {
		    if ($parts[2] ne "bacterium") {
			$strain = $parts[2]; # what if there is more
		    } elsif ($#parts >= 3) {
			$strain = $parts[3]; # what if there is more
		    }
		}
	    }
	}
    } else {
	if (($DA[0] ne "uncultured") && ($DA[0] ne "unidentified")) {
	    $genus = $DA[0];
	}
    }

    my $primerFound = 0; # if forwrd primer found, this is set to 1
                         # if reverse primer also found, then set to 2.
    my $primerFoundSeq = "";
    # Then process the primers and find the amplicons, if any
    if ($s =~ $forPrimer) { # look for forward primer
	$primerFound = 1;
	$primerFoundSeq = "F-";
	my $s2 = $'; #' # This is a neat trick to get the sequence 
	                # follwoing the matched primer. 
	                # Thus, removes everything before FOR primer
        if ($s2 =~ $revPrimer2) { # look for reverse primer
	    $primerFound = 2;
	    $primerFoundSeq = "FR";
	    my $amplicon = $`; # this is a neat trick to get the sequence 
                               # prior to the matched primer. 
	                       # Thus, removes everything after REV primer
	    # remember that the PCR fragment includes the 2 primers.
	    $amplicon = $amplicon.$revPrimer2;
	    $amplicon = $forPrimer.$amplicon;
	    my $l = length($amplicon);
	    # report length of PCR product and the name of the organism
	    print FH "$l\t$domain\t$genus\t$species\t$strain\t$gid\t$count\t$primerFoundSeq";
	    #if ($debug == 1) {
		print FH "\t$d\n";
	    #} else {
#		print FH "\n";
#	    }
	    # print lout "$count\t$primerFoundSeq\t$l\t$domain\t$genus\t$species\t$strain\t$gid\n";
        } elsif ($s2 =~ $revPrimer) { # look for reverse primer
	    $primerFound = 2;
	    $primerFoundSeq = "FF";
	    my $amplicon = $`; # this is a neat trick to get the sequence 
                               # prior to the matched primer. 
	                       # Thus, removes everything after REV primer
	    # remember that the PCR fragment includes the 2 primers.
	    $amplicon = $amplicon.$revPrimer;
	    $amplicon = $forPrimer.$amplicon;
	    my $l = length($amplicon);
	    # report length of PCR product and the name of the organism
	    print FH "$l\t$domain\t$genus\t$species\t$strain\t$gid\t$count\t$primerFoundSeq";
	    # print lout "$count\t$primerFoundSeq\t$l\t$domain\t$genus\t$species\t$strain\t$gid";
	    #if ($debug == 1) {
		print FH "\t$d\n";
	    #} else {
#		print FH "\n";
#	    }
        } else {
	    # print lout "$count\t$primerFoundSeq\t$l\t$domain\t$genus\t$species\t$strain\t$gid\n
	    #if ($debug == 1) {
		print FH "-\t$domain\t$genus\t$species\t$strain\t$gid\t$count\t$primerFoundSeq\n";
	    #} else {
		# print "\n";
	    #}
	}
    }
    print LFILE "\n";
}
   return;
}



