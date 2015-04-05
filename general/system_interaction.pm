#!/usr/bin/perl -w

package system_interaction;

use strict;

sub check_module_installed{
	eval("use $_[0]");
        if($@) {
		print STDERR "Module $_[0] is not installed\n";
                print STDERR "If you want to use it, please install $_[0]\n";
                print STDERR "Waiting for 10 seconds to continue\n";
                for(my $i = 1; $i < 10; $i++) {
                        print STDERR " . ";
                        sleep(1);
                }
                print STDERR "\n";
		return 1;
        }
	return 0;
}
sub check_multithreading_support{
	#Needs multithreadin support variable and user defined cores
        #First check if multithreading is not turned off by user   
	my $core_user = $_[1];     
        if($_[0] == 1) {
                return 1;
        }
        #Now check if multithreading is compiled into perl
        my $a = `perl -V | grep useithreads`;
        if(length($a) > 0) {
                #Threading is possible - check cores
                my $number = `nproc`;
                chomp $number;
                if($core_user == 0) {
                        return int($number/4);
                }
                if($core_user > $number) {
                        print STDERR "User specified more cores than exisiting!\n";
                        print STDERR "Use 1/4 of existing cores\n";
                        return int($number/4);
                } else {
                        return $core_user;
                }
        } else {
                return 1;
        }
        return 0;
}
1;
