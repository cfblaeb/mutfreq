# README #

### What is this repository for? ###

This project is for analysing mutations in amplicons


##### History and lessons learned
* History 
  * Project started out being named pparp after the first protein we studied. I then changed it to mutfreq and made it more generic. A lot of lessons regarding sequencing has been learned and a few of them I have written down.
  * Due to "the learning" there is probably some code left over that's no longer relevant. I try to keep the "Currently preferred method" up to date.

* Lessons
  * Remember that illumina reads are actually quite prone to SNPs. Some have reported error rates as high as 1 SNP every 500 basepairs.
  * Unless you really care about indels, you can save a lot of time by focusing on SNPs.
  * Remember that a 150bp read WT will, due to illumina errors, generate some high-ish count single snp mutants  

#### Currently preferred method:
  * Merge reads (e.g. with FLASH)
  * Filter using align.filter_flash_merged)
    * Remove anything with incorrect read length
    * Anything with a single base quality less than a threshold
  * Map (e.g. with bowtie2...maybe replace with pure python?)
  * Reads passing filter gets stored as their MD field in a counter object
  * Anything with a CIGAR string that doesnt match f"{merged_read_length}M" (but count them still for QA)
  
# NOTES TO SELF:
Looks like low quality scores are unavoidable....
Think about whether it makes sense to filter on quality or only on copy numbers....noise should be MUCH lower and you should be able to do that thing where you look at the jump in mean read count as you filter on higher and higher copy numbers.