
###
# Align the Mikkelsen MEF data with Bowtie and dump the resulting BAM files
# on S3
# Run Macs peak finding on them.

require 'catpaws'

#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, 'ami-20794c54'  #EC2 eu-west-1 64bit Lucid
set :instance_type, 'm1.large'
set :s3cfg, ENV['S3CFG'] #location of ubuntu s3cfg file
set :working_dir, '/mnt/work'
#ami-52794c26 32bit Lucid
#ami-505c6924 64bit Maverick
#ami-20794c54 64bit Lucid

set :nhosts, 1
set :group_name, 'mef_h3k4me3_mikkelsen_mef'

set :snap_id, `cat SNAPID`.chomp #ec2 eu-west-1 
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 50  #Needs to be the size of the snap plus enough space for alignments
set :ebs_zone, 'eu-west-1a'  #is where the ubuntu ami is
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'

#grab the latest snapshot ID for the raw Mikkelsen data (fastq)
task :get_snap_id, :roles=>:master do
  `curl http://github.com/cassj/mikkelsen07/raw/master/SNAPID > SNAPID `
end 


#make a new EBS volume from this snap 
#cap EBS:create

#and mount your EBS
#cap EBS:attach
#cap EBS:mount_xfs




# Bowtie v0.12.7
bowtie_link = 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie%2F0.12.7%2F&ts=1287377285&use_mirror=mesh'

#should really save this as an ami, but no time at the moment
desc "install bowtie"
task :install_bowtie, :roles => group_name do
  run "sudo apt-get update"
  run "sudo apt-get install -y zip unzip"
  run "cd #{working_dir} && wget -Obowtie.zip #{bowtie_link}" 
  run "cd #{working_dir} && unzip bowtie.zip"
  run "sudo cp #{working_dir}/bowtie*/bowtie* /usr/local/bin/"
end 
before 'install_bowtie', 'EC2:start'




# s3 setup
desc "install s3 client"
task :install_s3, :roles => group_name do
  sudo 'apt-get update'
  sudo 'apt-get install -y s3cmd'
end
before 'install_s3', 'EC2:start'

desc "upload the s3 config file"
task :s3_config, :roles => group_name do
  upload(s3cfg, '/home/ubuntu/.s3cfg')
end 
before 's3_config', 'EC2:start'


#get the current mouse genome (which I already have on S3).
task :fetch_genome, :roles => group_name do
  run "s3cmd get --force s3://bowtie.mm9/mm9.ebwt.zip #{working_dir}/bowtie-0.12.7/indexes/mm9.ebwt.zip"
  run "rm -Rf #{working_dir}/bowtie-0.12.7/indexes/chr*"
  run "cd  #{working_dir}/bowtie-0.12.7/indexes && unzip -o mm9.ebwt.zip"
end 
before "fetch_genome","EC2:start"


# run bowtie on the fastq file
# This is early Illumina data, I'm assuming qualities are solexa ASCII(QSolexa+64) scores. but they might not be. 
# they *might* be standard fastq 
task :run_bowtie, :roles => group_name do
  h3k4 = capture "ls #{mount_point}/MEF_H3K4me3_ChIPSeq/*fastq"
  h3k27 = capture "ls #{mount_point}/MEF_H3K27me3_ChIPSeq/*fastq"
  h3k4 = h3k4.split("\n")
  h3k27 = h3k27.split("\n")
  h3k4.each do |f|
    o = f.sub('.fastq', '.sam')
    run "#{working_dir}/bowtie*/bowtie --sam --al --best --solexa-quals  -q mm9 #{f} > #{o}"
  end
  h3k27.each do |f|
    o = f.sub('.fastq', '.sam')
    run "#{working_dir}/bowtie*/bowtie --sam --al --best --solexa-quals  -q mm9 #{f} > #{o}"
  end
end
before "run_bowtie", "EC2:start"



# fetch samtools from svn
desc "get samtools"
task :get_samtools, :roles => group_name do
  sudo "apt-get -y install subversion"
  run "cd #{working_dir} && svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools"
end
before "get_samtools", "EC2:start"


desc "build samtools"
task :build_samtools, :roles => group_name do
  sudo "apt-get -y install zlib1g-dev libncurses5-dev"
  run "cd #{working_dir}/samtools && make"
end
before "build_samtools", "EC2:start"


desc "install samtools"
task :install_samtools, :roles => group_name do
  sudo "cp #{working_dir}/samtools/samtools /usr/local/bin/samtools"
end
before "install_samtools", "EC2:start"



desc "make bam from sam"
task :to_bam, :roles => group_name do
  run "wget -O #{working_dir}/mm9_lengths  'http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths'"
  k4 = capture "ls #{mount_point}/MEF_H3K4me3_ChIPSeq"
  k4 = k4.split("\n").select{|f| f.match(/\.sam/)}
  k27 = capture "ls #{mount_point}/MEF_H3K27me3_ChIPSeq"
  k27 = k27.split("\n").select{|f| f.match(/\.sam/)}
  k4.each{|f| 
    f_out = f.sub('.sam', '.bam')
    run "samtools view -bt #{working_dir}/mm9_lengths -o #{mount_point}/MEF_H3K4me3_ChIPSeq/#{f_out} #{mount_point}/MEF_H3K4me3_ChIPSeq/#{f}"
  }
  k27.each{|f| 
    f_out = f.sub('.sam', '.bam')
    run "samtools view -bt #{working_dir}/mm9_lengths -o #{mount_point}/MEF_H3K27me3_ChIPSeq/#{f_out} #{mount_point}/MEF_H3K27me3_ChIPSeq/#{f}"
  }
end
before "to_bam", "EC2:start"


desc "sort bam"
task :sort_bam, :roles => group_name do
  k4 = capture "ls #{mount_point}/MEF_H3K4me3_ChIPSeq"
  k4 = k4.split("\n").select{|f| f.match(/\.bam/)}
  k27 = capture "ls #{mount_point}/MEF_H3K27me3_ChIPSeq"
  k27 = k27.split("\n").select{|f| f.match(/\.bam/)}
  k4.each{|f| 
    f_out = f.sub('.bam', '_sorted')
    run "samtools sort #{mount_point}/MEF_H3K4me3_ChIPSeq/#{f}  #{mount_point}/MEF_H3K4me3_ChIPSeq/#{f_out}"
  }
  k27.each{|f| 
    f_out = f.sub('.bam', '_sorted')
    run "samtools sort #{mount_point}/MEF_H3K27me3_ChIPSeq/#{f}  #{mount_point}/MEF_H3K27me3_ChIPSeq/#{f_out}"
  }
end
before "sort_bam", "EC2:start"

desc "remove duplicates"
task :rmdups, :roles => group_name do
  k4 = capture "ls #{mount_point}/MEF_H3K4me3_ChIPSeq"
  k4 = k4.split("\n").select{|f| f.match(/sorted\.bam/)}
  k27 = capture "ls #{mount_point}/MEF_H3K27me3_ChIPSeq"
  k27 = k27.split("\n").select{|f| f.match(/sorted\.bam/)}
  k4.each{|f| 
    f_out = f.sub('_sorted', '_sorted_nodups')
    run  "cd #{mount_point}/MEF_H3K4me3_ChIPSeq && samtools rmdup -s #{f} #{f_out}"
  }
  k27.each{|f| 
    f_out = f.sub('_sorted', '_sorted_nodups')
    run "cd #{mount_point}/MEF_H3K27me3_ChIPSeq && samtools rmdup -s #{f} #{f_out}"
  }
end
before "rmdups", "EC2:start"


#may need to do this with picard samtools doesn't manage
#header properly
desc "merge bam"
task :merge_bam, :roles => group_name do
  k4 = capture "ls #{mount_point}/MEF_H3K4me3_ChIPSeq"
  k4 = k4.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}.join(" ")
  k27 = capture "ls #{mount_point}/MEF_H3K27me3_ChIPSeq"
  k27 = k27.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}.join(" ")
  
  run "cd #{mount_point}/MEF_H3K4me3_ChIPSeq && samtools merge merged.bam #{k4}"
  run "cd #{mount_point}/MEF_H3K27me3_ChIPSeq && samtools merge merged.bam #{k27}"
end
before "merge_bam", "EC2:start"



desc "index bam files"
task :index, :roles => group_name do

  run "cd #{mount_point}/MEF_H3K4me3_ChIPSeq && samtools index merged.bam merged.bai"
  run "cd #{mount_point}/MEF_H3K27me3_ChIPSeq && samtools index merged.bam merged.bai"

end
before "index", "EC2:start"



desc "download bam files"
task :get_bam, :roles => group_name do
  `rm -Rf results/alignment/bowtie` #remove previous results
  `mkdir -p results/alignment/bowtie`
  dirs = capture "ls #{mount_point}"
  dirs = dirs.split("\n")
  dirs.each {|d|
    `mkdir -p results/alignment/bowtie/#{d}`
    files = capture "ls #{mount_point}/#{d}"
    files = files.split("\n").select{|f| f.match(/merged/)}
    files.each{|f|
      download( "#{mount_point}/#{d}/#{f}", "results/alignment/bowtie/#{d}/#{f}")
    }
  }
end
before "get_bam", 'EC2:start'

### Macs ?

macs_url ="http://liulab.dfci.harvard.edu/MACS/src/MACS-1.4.0beta.tar.gz"
macs_version = "MACS-1.4.0beta"

task :install_macs, :roles => group_name do
  sudo "apt-get install -y python"
  run "cd #{working_dir} && wget --http-user macs --http-passwd chipseq #{macs_url}"
  run "cd #{working_dir} && tar -xvzf #{macs_version}.tar.gz"
  run "cd #{working_dir}/#{macs_version} && sudo python setup.py install"
  sudo "ln -s /usr/local/bin/macs* /usr/local/bin/macs"
end
before "install_macs", 'EC2:start'

task :install_peaksplitter, :roles => group_name do
  url ='http://www.ebi.ac.uk/bertone/software/PeakSplitter_Cpp_1.0.tar.gz'
  filename = 'PeakSplitter_Cpp_1.0.tar.gz'
  bin = 'PeakSplitter_Cpp/PeakSplitter_Linux64/PeakSplitter'
  run "cd #{working_dir} && curl #{url} > #{filename}"
  run "cd #{working_dir} && tar -xvzf #{filename}"
  run "sudo cp #{working_dir}/#{bin} /usr/local/bin/PeakSplitter"
end 
before 'install_peaksplitter', 'EC2:start'

#you'll need to have done "install_r" and install_peak_splitter to do this
task :run_macs, :roles => group_name do
  
#  treatments = ["MEF_H3K27me3_ChIPSeq",
#               "MEF_H3K4me3_ChIPSeq"]
  treatments = ["MEF_H3K27me3_ChIPSeq"]

  control = "#{mount_point}/MEF_WCE_ChIPSeq/merged.bam"
  genome = 'mm'
  bws = [300]
  pvalues = [0.00001]

  #unsure what p values and bandwidths are appropriate, try a few.
  treatments.each{|t|
    treatment = "#{mount_point}/#{t}/merged.bam"
    bws.each {|bw|
      pvalues.each { |pvalue|
        
        name = "#{t}"
        dir = "#{mount_point}/macs_#{bw}_#{pvalue}_#{t}"
        run "rm -Rf #{dir}"
        run "mkdir #{dir}"
        
        macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{name} --format BAM --gsize #{genome} --bw #{bw} --pvalue #{pvalue}"
        run "cd #{dir} && #{macs_cmd}"
        
        dir = "#{mount_point}/macs_#{bw}_#{pvalue}_#{t}_subpeaks"
        run "rm -Rf #{dir}"
        run "mkdir #{dir}"
        
        # With SubPeak finding
        # this will take a lot longer as you have to save the wig file 
        macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{group_name} --format BAM --gsize #{genome} --call-subpeaks  --bw #{bw} --pvalue #{pvalue} --wig"
        run "cd #{dir} && #{macs_cmd}"
        
      }
    }
  }
  
end
before 'run_macs', 'EC2:start'



#pack up the runs and downloads them to the server (without the wig files)
task :pack_macs, :roles => group_name do
  macs_dirs = capture "ls #{mount_point}"
  macs_dirs = macs_dirs.split("\n").select {|f| f.match(/.*macs.*/)}
  macs_dirs.each{|d|
    run "cd #{mount_point} &&  tar --exclude *_wiggle* -cvzf #{d}.tgz #{d}"
  }
  
end
before 'pack_macs','EC2:start' 

task :get_macs, :roles => group_name do
  macs_files = capture "ls #{mount_point}"
  macs_files = macs_files.split("\n").select {|f| f.match(/.*macs.*\.tgz/)}
  res_dir = 'results/alignment/bowtie/peakfinding/macs'
  `rm -Rf #{res_dir}`
  `mkdir -p #{res_dir}`
  macs_files.each{|f| 
    download("#{mount_point}/#{f}", "#{res_dir}/#{f}") 
    `cd #{res_dir} && tar -xvzf #{f}`
  }

end
before 'get_macs', 'EC2:start'





#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




