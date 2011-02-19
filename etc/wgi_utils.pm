package wgi_utils;

sub which_bkgd
{
    @lines_in=@_;

#   We want the last model read if there is more than one
    my @lines = reverse sort(@lines_in);

# Determine background model.

    $modeltype = "missing";

#Newer way with using 'get_modelfg_3d'
    $searchstring = "Successfully obtained: ";
    foreach (@lines)
    {
       if (/$searchstring(.*)/)
       {
           ($modeltime,$modeltype) = split(" ",$1);
           last;
       }
    }

    if ($modeltype eq "missing"){
        print "Can't determine model background.\n";
        return;
    }

    if ($modeltype eq 'fua' || $modeltype eq 'fsf') {
        $searchstring = "Reading 3d";
        $searchstring = "Reading 2d" unless ($modeltype eq 'fua');
        foreach (@lines)
        {
           if(/$searchstring(.*)/)
           {
               @string_elems=split(" ",$1);
               $bkgdpath = @string_elems[0];
               @string_elems = split("/",$bkgdpath);
               $modelid = @string_elems[$#string_elems];
               last;
           }
        }
    }

    print "BACKGROUND FIELDS:  ";

    if ($modeltype eq 'lga' || $modeltype eq 'lgb')
    {

        $bgmodelfile = $main::filename;
        $bgmodelfile =~ s/$main::logname/lga/;
        $modelid = "???";
        if(open(BGMODEL,$bgmodelfile))
        {
           foreach(<BGMODEL>)
           {

              if(/cmodel (\w*)/) {$modelid = $1;}
              if (/cmodel: (.*)/)           {(@modelid) = split(" ", $1); $modelid = @modelid[0];}
              if (/Reading - (.*)/)         {$pathname = $1; last; close BGMODEL; }
              if (/reading cdfname\: (.*)/) {$pathname = $1; last; close BGMODEL; }
           }
           $basename = $1 if ($pathname =~ /([^\/]*)$/);
        }

        if($basename =~ /(\d\d\d\d\d)(\d\d\d\d)(\d\d)(\d\d)/)
        {
           $runtime = $2/100;
           $fcsthr = $4;
#   print "Using $fcsthr hr fcst from $modelid started at $runtime UTC \n";
        }
        elsif($basename =~ /(\d\d\d\d)(\d\d)(\d\d)\_(\d\d)(\d\d)/)
        {
              $runtime = $4;
#   print "Using grids from $modelid started at $runtime UTC \n";
        }

    }

    if($modeltime =~ /(\d\d)(\d\d\d)(\d\d\d\d)(\d\d\d\d)/ )
    {
       $runtime = $3;
       $fcsthr = $4;
       $lengthmdlid=length $modelid;
#      print "Using $fcsthr fcst from $modelid model started at $runtime UTC \n";
#      print "$modeltype $modelid $1$2$3$4 \n";
    }

    return ($modelid,$modeltype,$runtime,$fcsthr,$1,$2);

}
1;
#
#--------------------------------------------------------------
#
sub get_log_filename
{
    my ($dataroot,$log_name,$hour)=@_;
    if(!defined $dataroot){$dataroot = "./..";}
    my @fnames = <$dataroot/log/$log_name.log*>;
    my $nf = @fnames;
    if ($nf == 0) {die "No log files found.\n";}

    open(SYSDAT,"$dataroot/time/systime.dat");
    my @systime_info=<SYSDAT>;
    close SYSDAT;
    my $yyjjjhhmm = $1 if($systime_info[1] =~ /(\d\d\d\d\d\d\d\d\d)/);
    my ($file,$fname_yyjjjhhmm,$fname);
    my $age_of_newest = 100;
    foreach $file (@fnames){
            $age = -M $file;
            if ($age < $age_of_newest){
                $fname = $file;
                $fname_yyjjjhhmm = $1 if($fname =~ /(\d\d\d\d\d\d\d\d\d)/);
                $age_of_newest = $age;
                
            }
    }

# This is a bit more kosher... this doesn't work with log file names having ".log.yyjjjhhmm"
    if (defined $hour) { 
           $fname = $dataroot."/log/".$log_name.".log.".$hour;
    }else{
           if($fname_yyjjjhhmm != $yyjjjhhmm){
              $fname = $dataroot."/log/".$log_name.".log.".$yyjjjhhmm;
           }
    }
    return $fname;
}
1;
