# getParentOfOrigin
The script to detect parent-of-origin for denovo mutations by CRAM 


### Usage:
```
perl getParentOfOrigin_hg38_cram.pl -i ./Example_input.txt -c XXX.cram 
```
- `-i` pairs of denovo mutation and informative variant in the same individual (Please add a path)
- `-c` individual cram. The index file (`XXX.cram.crai`) should be stored in the same folder


### Input file format
#### Columns are listed below (no need to add header line)
```
distance	chr_dn	pos_dn	ref_dn	var_dn	sample	SSCID	fam	type	chr_info	pos_info	ref_info	var_info	info_origin	info_gt_code
```
- distance: distance between variant pairs
- info_origin: Origin of informative variant: `Fa` (paternal) or `Mo`(maternal)
- (optional) info_gt_code: Genotype code for informative variant (Order: father/mother/child; For example: 101):
```
0/0 - 0
0/1 - 1
1/1 - 2
```
- Note: Key columns are: chr/pos/ref/var for both variants, as well as info_origin. The rest columns should be present to run the script (but they can be deleted by tiny modifications)

### Output file format 
- The result line will be generated when:
    - Observed at least one overlap read/read pair
    - The proportion of the variant reads (PVR) for denovo mutation (heterozygous) > 10% 
    - The PVR for informative variant (heterozygous) between 10% - 95%
    - Please note these are very simple filters for heterozygous variants (other HQ filters were done before)
- Result summary for each denovo mutation: `Checking_results.` + input_prefix (`Example_input`)  + `.txt` (no header line)
#### New columns will be added to the summary file (output file no header line):
```
result	confidence_level	(input columns) dn_total_count	dn_support_count	percent_dn	info_total_count	info_support_count	percent_info	overlap_count	consist_count	percent_consist	incosist_count	percent_inconsist
```
- result: 
    - Fa/Mo: decision based only the consistency outcome on overlapped reads/read pair
    - Unclear: if the percent of consistent/inconsistent < 90%
- confidence_level
    - Confident: If overlapped reads/read pair > 5
    - Weak: If overlapped reads/read pair <= 5



### INSTALLATION 
1. Install htslib
(reference: https://github.com/samtools/htslib)
```
git clone -b master --depth=1 https://github.com/samtools/htslib.git
cd htslib
git submodule update --init --recursive
sudo make install
sudo ldconfig
cd ..
```
2. Install Bio-DB-HTS
(reference: https://github.com/Ensembl/Bio-DB-HTS)
```
git clone -b master https://github.com/Ensembl/Bio-DB-HTS.git
cd Bio-DB-HTS
perl Build.PL
sudo ./Build install
cd t
for f in $(ls *.t) ;
    do
        perl $f
    done
cd ..
```
- Test installation:
should see no error when running the following command
```
perl test_installation.pl
```
### TROUBLESHOOTING (optional)

- (To save running time) Note for samtools (CRAMs): The reference files (`.fa` and `.fai`) should be located in the same path as showed in the header
- For example: `/data/NYGC/Resources/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa`
```
[W::cram_populate_ref] Creating reference cache directory /root/.cache/hts-ref
This may become large; see the samtools(1) manual page REF_CACHE discussion
```
