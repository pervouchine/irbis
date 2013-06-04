#!/usr/bin/perl

$pid = fork();

if( $pid == 0 ){
    exec(@ARGV);
    exit(0);
}

$retval = time();
$memory = 0;
while(1) {
    ($id, $usage) = split /\n/, `ps -p $pid -o size`;
    last unless($usage>0);
    $memory = $usage if($memory<$usage);
}

open FILE, ">>memtest.log";
print FILE "$retval\t$memory\t", time() - $retval, "\t",join(" ", @ARGV),"\n";
close FILE;
