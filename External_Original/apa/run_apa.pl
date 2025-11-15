$f0=$ARGV[0];
open(F,$f0);
while(<F>){
    ($batch,$sample,$bam)=split;
    $id=$batch."_".$sample;
    open(SL,"> ${id}.slum");
    print SL  "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=4\n#SBATCH --mem=50G\n#SBATCH --time=1-00:01:00\n#SBATCH --mail-user=zhishan.chen\@vanderbilt.edu\n#SBATCH --mail-type=ALL\n#SBATCH -o  $id.log\n";
    print SL "                                \n";
    print SL "ml GCC/8.2.0 BEDTools/2.28.0 \n";
    print SL "bedtools genomecov -ibam $bam -bga -split -trackline \> wig/${sample}.wig \n";
    close(SL);
    print "${id}:";
    sleep 10;
    system "sbatch ${id}.slum";
}

