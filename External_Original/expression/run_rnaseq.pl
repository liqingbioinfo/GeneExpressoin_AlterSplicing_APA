$f0=$ARGV[0];
open(F,$f0);
while(<F>){
    ($batch,$id,$sample)=split;
    $id2=$batch."_".$id;
    open(SL,"> ${id2}.slum");
    print SL  "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=50G\n#SBATCH --time=1-01:00:00\n#SBATCH --mail-user=zhishan.chen\@vanderbilt.edu\n#SBATCH --mail-type=ALL\n#SBATCH -o  $id2.log\n";
    print SL "                                \n";
    print SL "bash rnaseq.sh $batch $sample \n" ;
    close(SL);
    print "${id2}:";
    sleep 1;
    system "sbatch ${id2}.slum";
}

