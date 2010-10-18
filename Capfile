
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

#I am here

desc "sort bam"
task :sort_bam, :roles => group_name do
  run "samtools sort #{working_dir}/export.bam #{working_dir}/export_sorted"
end
before "sort_bam", "EC2:start"

#merge bam files here!

desc "remove duplicates"
task :rmdups, :roles => :chipseq do
  run "samtools rmdup -s #{working_dir}/export_sorted.bam #{working_dir}/export_nodups.bam"
end
before "rmdups", "EC2:start"

desc "index bam files"
task :index, :roles => :chipseq do
  run "samtools index #{working_dir}/export_nodups.bam #{working_dir}/export_nodups.bai"
end
before "index", "EC2:start"



desc "upload to s3"
task "to_s3", :roles => :chipseq do

  servers = find_servers_for_task(current_task)
  puts servers
  it = 0..(servers.length - 1)
  it.each do |i|
    host = servers[i]
    bucket = bucket_names[i]
    object = object_names[i]
    bucket = "bam."+bucket
    run("s3cmd mb s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export.bam s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export_sorted.bam s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export_nodups.bam s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export_nodups.bai s3://#{bucket}", :hosts => host)
  end

end
before "to_s3", "EC2:start"





#if you want to keep the results

#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




