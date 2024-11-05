########################################################################################
线粒体基因组组装

0. 获得测序数据后，需要
1) 查看测序报告，用md5sum或其他校验工具检查数据完整性
2) fastqc/multiqc评估测序数据质量
3) 如果测序报告中没有污染评估，需要自己评估，参考https://blog.csdn.net/qq_42962326/article/details/105081327；但这篇文章中涉及blast '-max_target_seqs 1'这一参数的问题，需要修改；修改办法是没有用'-max_target_seqs 1'这个参数，blast生成的xml格式的，用survey-tiqu.py脚本把每一个匹配结果的第一条选出来，再进行后面的步骤；这一部分晓琦做过很多，可以问她


1. install mitoz version 2.3
# the installation may be out-of-date, check https://github.com/linzhi2013/MitoZ/wiki/Installation
# download this file first: https://github.com/linzhi2013/MitoZ/blob/master/version_2.3/mitozEnv.yaml
# create the environment mitozEnv for mitoz, and install dependancies for mitoz in mitozEnv.yaml
conda env create -n mitozEnv -f mitozEnv.yaml
conda activate mitozEnv
conda install -c bioconda mitoz


2. install NCBI taxonomy database
python3
>>> from ete3 import NCBITaxa
>>> ncbi = NCBITaxa()
>>> ncbi.update_taxonomy_database()
# ERROR: sqlite3.IntegrityError: UNIQUE constraint failed: synonym.spname, synonym.taxid
# SOLUTION: https://github.com/linzhi2013/MitoZ/issues/72#issuecomment-666215301
# 链接:https://pan.baidu.com/s/1Sa8icV5TVNhfKR4PHWmgow 密码:asuk
# Download the file 'database.etetoolkit.tar.bz2' into your linux machine, place it into your home directory, i.e. '~/database.etetoolkit.tar.bz2'.
cd 
tar -xjf database.etetoolkit.tar.bz2
# will get a .etetoolkit directory
# test if the database works fine
python3
>>> from ete3 import NCBITaxa
>>> a = NCBITaxa()
>>> a.get_name_translator(["Arthropoda"])
{'Arthropoda': [6656]}


3. extract a subset of reads
# 3Gb/end is enough
1) 'extractfq' # extractfq 提取前i Gb序列
pip install extractfq
extractfq -fq1 /data/backup/2021_Tbrevicauda_Tstewarti_Tstenura/clean/Tlepto-01F/Tlepto_01F_HGCM5DSX2_L4_1.clean.fq.gz -fq2 /data/backup/2021_Tbrevicauda_Tstewarti_Tstenura/clean/Tlepto-01F/Tlepto_01F_HGCM5DSX2_L4_2.clean.fq.gz -outfq1 Tlepto_01F_1_extractfq.clean.fq.gz -outfq2 Tlepto_01F_2_extractfq.clean.fq.gz -size_required 2 -gz
# 看起来2 Gb有时不够用，需要3 Gb或更多，我记得-size_required设置的是双端数据一共几G

2) 'seqtk' # seqtk sample 是随机subsample序列
seqtk sample -s100 ${FILE1} 20000000 > ${OUTPUT1}
seqtk sample -s100 ${FILE2} 20000000 > ${OUTPUT2}


4. run MitoZ 
#参考说明https://gitee.com/CHANyp/MitoZ#63-directory-structure
#生成的.megahit.result目录中的.megahit.mitogenome.fa文件是线粒体基因组文件，将用于后续分析
#大约1.5至3G碱基对（bp）足以进行线粒体基因组组装
source activate mitozEnv
python3 /home/wangyf/MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha/MitoZ.py all \
--genetic_code 2 --clade Chordata \
--outprefix Esox_lucius_BF301-02R0001_good_mitoz \
--thread_number 20 \
--fastq1 /home/wangyf/wyf/bmkdata/data/Unknown_BF301-02R0001_good_1.fq.gz \
--fastq2 /home/wangyf/wyf/bmkdata/data/Unknown_BF301-02R0001_good_2.fq.gz \
--fastq_read_length 150 --insert_size 350 \
--run_mode 2 \
--filter_taxa_method 1 \
--requiring_taxa 'Chordata'
--species_name Esox_lucius


5. reorder mitogenome
# Mitogenome_reorder.py将基因组重新排列，该文件在/home/wangyf/MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha/useful_scripts/里面
# 需要一个reference mitogenome fasta 文件且开头>中必须含有topology=circular|linear等字样
awk 'BEGIN{print ">Trosa ref topology=circular"} NR>1{print}' Trosa_ref_mitogenome.fasta > Trosa_ref_mitogenome_circular.fasta
python3 /home/wangyf/MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha/useful_scripts/Mitogenome_reorder.py -f Tstoli_01F_mito.result/work71.mitogenome.fa -r ../Trosa_ref_mitogenome_circular.fasta && mv work71.mitogenome.fa.reorder Tstoli_01F_mitoz_reorder.fa


6. 也可使用GetOrganelle v1.7.5 有参组装线粒体序列，*.path_sequence.fasta是组装结果文件
#输出目录需要是一个不存在的目录
#outfile——*.path_sequence.fasta, each fasta file represents one type of genome structure；如果生成的基因组是完整的（名称为*.fasta），则可以删除除此之外的文件。

conda install -c bioconda getorganelle
get_organelle_config.py --add animal_mt
get_organelle_from_reads.py -1 Unknown_BF301-03R0010_good_1.fq.gz -2 Unknown_BF301-03R0010_good_2.fq.gz -R 10 -k 21,45,65,85,105 -F animal_mt -o jianqiju/jianqiju-out


7. MITOS2网页版预测线粒体蛋白编码序列 # MITOS2也可以预测非编码序列
# 序列较多可以合成一个dataset一起提交
# 序列数目特别多的话，MITOS2可能不太方便，mitoz也可以在组装后预测，但是我用的mitoz anaconda版本预测出来结果会很奇怪，只能预测出来2、3个蛋白编码基因，你可以再试试你的版本，也可以找找其他更方便的注释软件


###################################################################################################################################################
线粒体基因组注释，使用mitofish网页


###################################################################################################################################################
## 建树

#需要软件trimal、macse、translatorx、samtools、seqkit、iqtree
#如果有的基因缺失，需要将缺失的序列用NNN补上，否则串联时该个体的其他序列也串联不起来，造成个体数据缺失
# 对每个基因单独进行比对和过滤
# 因为translatorx运行完毕后会将序列的原顺序打乱，因此需要提前将序列命名为物种名，以便后续用seqkit concat同一物种不同基因的序列串联在一起
# translatorx和trimal安装在 conda trans 环境下，且需要translatorx这个文件的路径，复制了一份在 home/wangyf/translatorx

cd /home/wangyf/wyf/bmkdata/liyu-alignment/83gene

source activate trans
cat mito_gene_list | while read gene; do \
/home/wangyf/translatorx -i ${gene}.fa -o ${gene}.transx -p F -c 2; \ # 按照编码的氨基酸进行序列比对
java -jar /home/wangyf/macse_v2.06.jar -prog exportAlignment -align ${gene}.transx.nt_ali.fasta -codonForInternalStop NNN -codonForFinalStop --- -gc_def 02; \ # 对终止密码子进行标注，后面是对生成的NT文件建立索引
trimal -in ${gene}.transx.nt_ali_NT.fasta -out ${gene}.transx.macse.trimal -automated1; \ # 使用trimal进行比对的trimming。先肉眼看一眼，如果比对的不错，这一步可跳过。
/home/wangyf/translatorx -i ${gene}.transx.macse.trimal -o ${gene}.transx.macse.trimal.transx.fasta -p F -c 2;\ # 修剪后再次进行序列比对
java -jar /home/wangyf/macse_v2.06.jar -prog exportAlignment -align ${gene}.transx.macse.trimal.transx.fasta.nt_ali.fasta -codonForInternalStop NNN -codonForFinalStop --- -gc_def 02;
done
# 注意对生成的比对文件进行检查，比如比对是否有异常，长度是否正常
#此处需要修改文件里的ID名,将每个基因文件中的每个样本名改成一样的，后面才能串联到一起
sed -i 's/[_-]ATPase8/ /g' atp8.152samples.transx.nt_ali_NT.fasta
# 对每个基因重复操作，或想办法改成一个循环

#####
# 对生成的fasta文件建立索引

conda activate samtools
cat mito_gene_list| while read gene; do samtools faidx ${gene}.transx.macse.trimal.transx.fasta.nt_ali_NT.fasta; done

# mito_gene_list文件列出线粒体基因名
cat mito_gene_list| while read gene; do ls ${gene}.transx.macse.trimal.transx.fasta.nt_ali_NT.fasta; done > infile.list

# 将各个基因的alignment串联到一起
seqkit concat --infile-list infile.list --out-file concat.fasta

# 分partition构建线粒体蛋白编码基因系统发生树
# 制作一个partition file，参考如下命令，将每个氨基酸的第1个碱基、第2个碱基、第3个碱基分成不同的partition
cat mito_gene_list | while read gene; do echo ${gene}; head -n 1 ${gene}.transx.macse.trimal.transx.fasta.nt_ali_NT.fasta.fai | awk '{print $2}'; done | paste - - | awk 'BEGIN{tot=0; line=0} {print "DNA, part"line*3+1" = "tot+1"-"tot+$2"\\3\n""DNA, part"line*3+2" = "tot+2"-"tot+$2"\\3\n""DNA, part"line*3+3" = "tot+3"-"tot+$2"\\3";tot=tot+$2;line+=1}' > partition.file
# 建树
source activate iqtree
iqtree -s concat.fasta --seed 6574744 -T 5 -m MFP+MERGE -B 1000 --prefix cdstest -p partition.file

#得到的结果文件中.contree 和.treefile 文件下载下来，用FigTree软件/iTOL网站看

conda activate base
cat mito_gene_list| while read gene; do ls ${gene}.transx.macse.trimal.transx.fasta.nt_ali_NT.fasta; done > infile.list
seqkit concat --infile-list infile.list --out-file concat.fasta
cat mito_gene_list | while read gene; do echo ${gene}; head -n 1 ${gene}.transx.macse.trimal.transx.fasta.nt_ali_NT.fasta.fai | awk '{print $2}'; done | paste - - | awk 'BEGIN{tot=0; line=0} {print "DNA, part"line*3+1" = "tot+1"-"tot+$2"\\3\n""DNA, part"line*3+2" = "tot+2"-"tot+$2"\\3\n""DNA, part"line*3+3" = "tot+3"-"tot+$2"\\3";tot=tot+$2;line+=1}' > partition.file
cat mito_gene_list | while read gene; do echo ${gene}; head -n 1 ${gene}.transx.nt_ali_NT.fasta.fai | awk '{print $2}'; done | paste - - | awk 'BEGIN{tot=0; line=0} {print "DNA, part"line*3+1" = "tot+1"-"tot+$2"\\3\n""DNA, part"line*3+2" = "tot+2"-"tot+$2"\\3\n""DNA, part"line*3+3" = "tot+3"-"tot+$2"\\3";tot=tot+$2;line+=1}' > partition.file

source activate iqtree
iqtree -s concat.fasta --seed 6574744 -T 5 -m MFP+MERGE -B 1000 --prefix rna -p partition.file
