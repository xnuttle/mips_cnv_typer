//Xander Nuttle
//detail_mip_targets_v3.c
//Call: ./detail_mip_targets_v3 mip.armlocs (int)number_of_gene_families (int)number_of_contigs_in_1st_gene_family <(int)number_of_contigs_in_Nth_gene_family> master_sequence_for_1st_family.fasta <master_sequence_for_Nth_family.fasta> 1st_family_1st_contig.fasta 1st_family_2nd_contig.fasta <1st_family_Nth_contig.fasta> 1st_family.suns.fixed 1st_family_exons.bed <Nth_family_first_contig.fasta> <Nth_family_Nth_contig.fasta> <Nth_family.suns.fixed> <Nth_family_exons.bed> output_file_base_name
//
//This program differs from detail_mip_targets.c in that it uses the newer 5-column format for the SUNs file input (see below) and that it does not
//consider SUN score information, instead assuming the SUNs input were not actually scored using high-coverage genomes. It also does a smarter search for
//the MIP targets within contig sequences, considering both MIP arms and the proper distance between them instead of just a single MIP arm.
//
//This program differs from detail_mip_targets_v2.c in that MIP specificity in this program's output is encoded by a string of letters rather than 0s and 1s,
//where the string of letters is generated based on different distinct MIP target sequences. For example, if there are 4 paralogs and the first 2 paralogs have
//identical sequence over a MIP target, the last 2 paralogs have identical sequence over the same MIP target, but the sequences of the first and last paralog differ,
//the specificty for that MIP would be defined as "AABB", even though no paralogs have a distinct target sequence (i.e. the specificity using the v2 program would be
//encoded as "0000" because no target sequences are distinct from all others); this new description of MIP paralog-specificity is particularly useful for analyses
//considering multiple alleles of the same paralog, where only rarely (i.e. where there is a SNP) would allelic target sequences differ from each other, but where
//often all allelic sequences are distinct from all paralogous sequences.
//
//Newer SUN file format example:
//CH17_chr16_hydin        154911  380     0.000000        SUN_is_masked
//CH17_chr16_hydin        155382  851     11.800000       no_masking_near_SUN
//CH17_chr16_hydin        155567  1036    10.733333       no_masking_near_SUN
//CH17_chr16_hydin        155982  1451    11.566667       no_masking_near_SUN
//CH17_chr16_hydin        156053  1522    6.766667        SUN_near_masked_bases
//CH17_chr16_hydin        156066  1535    3.666667        SUN_near_masked_bases
//CH17_chr16_hydin        156934  2403    11.466667       no_masking_near_SUN
//CH17_chr16_hydin        157056  2525    0.000000        SUN_is_masked
//CH17_chr16_hydin        157111  2580    0.000000        SUN_is_masked
//
//This program generates a file containing information on MIP targets and prints out information on any targeting conflicts.
//
//Example call: ./detail_mip_targets SRGAP2_RH_mip.armlocs 2 4 2 SRGAP2_1q32.fasta RH_master.fasta ~/SUNK_analysis_SRGAP2_chr1q21_final/q32_contig_final.fasta ~/SUNK_analysis_SRGAP2_chr1q21_final/q21_contig_new_oct2011.fasta ~/SUNK_analysis_SRGAP2_chr1q21_final/p12_contig_new_oct_2011.fasta ~/SUNK_analysis_SRGAP2_chr1q21_final/CH17-266P3.fasta SRGAP2.suns.fixed SRGAP2_1q32_exons.bed ../RH/RHD_contig.fasta ../RH/RHCE_contig.fasta RH.suns.fixed RH_exons.bed SRGAP2_RH
//
//In order to analyze mapped sequence data from a MIP experiment, it is first necessary to detail MIP targets and store important information about them
//in a format that allows for straightforward input into an analysis program. This program performs these tasks for an arbitrary set of MIPs. Inputs are:
// -an armlocs file containing MIP targeted gene family(ies), and extension and ligation arm start and end coordinates (base 1) with respect to the
// targeted gene family's master sequence (this file can be easily created from information in the MIP design Excel file); MIPs in this file should be
// sorted by the gene family they target
// -a series of at least two integers, the first being the number of gene families targeted and the rest corresponding to numbers of contig sequences
// in each targeted gene family 
// -fasta files containing each gene family's master sequence (if there are multiple targeted gene families, the order of fasta files should match the
// order of integers specifying the number of contigs in each gene family)
// -for each gene family targeted, fasta files containing contig sequences, a fixed SUNs file, and a bed file specifying exon coordinates (base 1)
// in the master sequence
// -the base name of your desired output file (to which ".miptargets" will be appended to give you your final output file)
//
//Example line from input file "mip.armlocs" followed by an explanation of each field in brakets (each MIP gets its own line):
//SRGAP2_1q32 4543  4561  4410  4430
//[target_master_sequence extension_probe_start extension_probe_stop ligation_probe_start ligation_prob_stop]
//
//This program output a file having 9 fields containing information on each MIP target:
//  1-name of master sequence targeted
//  2-start and end coordinates (start,end) of targeted sequence plus arms with respect to master sequence
//  3-string of contig names and start coordinates of targeted sequence plus arms with respect to each contig (name:start;)
//  4-MIP type (S=SUN-targeting, E=exon-targeting, B=both)
//  5-string of locations of paralog-specificity-conferring variation (all such variation seen in aligned sequences) with respect to targeted sequence plus arms
//  	(loc1,loc2,loc3,); this will be used to flag such bases for base call quality assessment; the range of locations is 1-152 becaue targeted sequence plus arms is 152 bp
//  6-same as field 5, except only containing locations of fixed SUNs with respect to targeted sequence plus arms; this will be used to flag such bases for 
//  	assessment of mismatches at those locations
//  7-string of ones and zeros having length equal to the number of contigs corresponding to the targeted master sequence, with a "1" indicating the
//    corresponding targeted contig sequence is distinguishable from other contig sequences and a "0" indicating the opposite (1010); this example means paralogs
//    1 and 3 have unique sequence over the 152 targeted bases but paralogs 2 and 4 are indistinguishable
//  8-mip orientation with respect to master sequence (+ or -)
//	9-length of the first MIP arm in bp (first with respect to the reference contig coordinates)
//
//Example line from output file "SRGAP2_RH.miptargets"
//SRGAP2_1q32	4893,5044	1q32.1_final:123878;q21_contig_new_oct_2011:145255;p12_contig_new_oct_2011:275550;SRGAP2D_contig:59940;	S	34,41,103,	41,103,	1010	-	22
//
//Finally, this program also prints to standard output any targeting conflicts detected (instances of more than 1 MIP targeting the same strand of the same genomic region).
//Such situations are possible if MIP design for exonic regions and for SUN regions was not integrated to avoid such conflicts. Interpretaion of any data
//from MIPs having targeting conflicts should take the possibility of poor capture performance (and resulting poor data) due to a targeting conflict into account.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc,char*argv[])
{
	//get information from command line on number of gene families
	int num_gene_families;
	num_gene_families=strtol(*(argv+2),NULL,10);

	//get information from command line on number of contigs in each gene family and calculate maximum number of contigs in a family
	int num_contigs_per_family[num_gene_families];
	int i;
	int max_num_contigs=0;
	for(i=0;i<num_gene_families;i++)
	{
		num_contigs_per_family[i]=strtol(*(argv+3+i),NULL,10);
		if(num_contigs_per_family[i]>max_num_contigs)
		{
			max_num_contigs=num_contigs_per_family[i];
		}
	}

	//set up mip target structure
	struct mip_target
	{
		char target_sequence[51];
		char master_loc[31];
		char contig_start_locs[61*max_num_contigs];
		char mip_type;
		char mip_orientation;
		char variant_locs[601];
		char fixed_sun_locs[301];
		char specificity[max_num_contigs+1];
		int first_arm_length;
	};	

	//get names of master sequences
	char master_sequence_names[num_gene_families][51];
	FILE*master_sequence_files[num_gene_families];
	char ch;
	int j;
	for(i=0;i<num_gene_families;i++)
	{
		j=0;
		master_sequence_files[i]=fopen(*(argv+3+num_gene_families+i),"r");
		while(((ch=getc(master_sequence_files[i]))!='\n')&&(ch!='\r'))
		{
			if(ch!='>')
			{
				master_sequence_names[i][j]=ch;
				j++;
			}
		}
		master_sequence_names[i][j]='\0';
	}	


	//read in information about mip arm location and calculate total number of mips and number of mips targeting each family
	int num_mips=0;
	int num_mips_per_family[num_gene_families];
	for(i=0;i<num_gene_families;i++)
	{
		num_mips_per_family[i]=0;
	}
	FILE*mip_arm_locations_file;
	fpos_t pos;
	char master_sequence_target[51];
	master_sequence_target[50]='\0';
	long ext_start,ext_end,lig_start,lig_end;
	mip_arm_locations_file=fopen(*(argv+1),"r");
	fgetpos(mip_arm_locations_file,&pos);
	while(fscanf(mip_arm_locations_file,"%s %ld %ld %ld %ld",master_sequence_target,&ext_start,&ext_end,&lig_start,&lig_end)==5)
	{
		num_mips++;
		for(i=0;i<num_gene_families;i++)
		{
			if(strncmp(master_sequence_target,master_sequence_names[i],strlen(master_sequence_target))==0)
			{
				num_mips_per_family[i]++;
			}
		}
	}
	fsetpos(mip_arm_locations_file,&pos);

	//allocate memory to store inportant mip target information
	struct mip_target*mip_targets;
	mip_targets=(struct mip_target*)malloc(num_mips*sizeof(struct mip_target));
	if(mip_targets==NULL)
	{
		printf("Memory allocation failed!\n");
		return 1;
	}


	//for each gene family, read in information on master sequence, contig sequences, good SUNs, and exons
	//process mip arm locations file simultaneously and determine output fields for each mip target
	int num_prev_gene_family_files;
	int k=0,m=0,n=0;
	for(i=0;i<num_gene_families;i++)
	{
		//input master sequence
		char*master_sequence;
		master_sequence=(char*)malloc(1000001*sizeof(char)); //1 Mb size limit on master sequence
		if(master_sequence==NULL)
		{
			printf("Memory allocation failed!\n");
			return 1;
		}
		k=0;
		while((ch=getc(master_sequence_files[i]))!=EOF)
		{
			if((!(isspace(ch)))&&(isalpha(ch)))
			{
				master_sequence[k]=ch;
				k++;
			}
		}
		master_sequence[k]='\0';
		fclose(master_sequence_files[i]);
		
		//input contig sequences
		num_prev_gene_family_files=0;
		for(j=0;j<i;j++)
		{
			num_prev_gene_family_files+=num_contigs_per_family[j]+2; //number of contig sequence files + fixed SUNs file + exons file
		}
		char contig_sequence_names[num_contigs_per_family[i]][51];
		FILE*contig_sequence_files[num_contigs_per_family[i]];
		char*contig_sequences[num_contigs_per_family[i]];
		for(j=0;j<num_contigs_per_family[i];j++)
		{
			k=0;
			contig_sequence_files[j]=fopen(*(argv+3+2*num_gene_families+num_prev_gene_family_files+j),"r");
			while(((ch=getc(contig_sequence_files[j]))!='\n')&&(ch!='\r'))
    	{
      	if(ch!='>')
      	{
        	contig_sequence_names[j][k]=ch;
        	k++;
      	}
    	}
			contig_sequence_names[j][k]='\0';
			contig_sequences[j]=(char*)malloc(1000001*sizeof(char));
			if(contig_sequences[j]==NULL)
			{
				printf("Memory allocation failed!\n");
				return 1;
			}
			k=0;
			while((ch=getc(contig_sequence_files[j]))!=EOF)
			{
				if((!(isspace(ch)))&&(isalpha(ch)))
				{
					contig_sequences[j][k]=ch;
					k++;
				}
			}
			contig_sequences[j][k]='\0';
			fclose(contig_sequence_files[j]);
		}

		//input suns !!! (in this version [detail_mip_targets_v2.c] information on all SUNs is input, not only SUNs considered to be fixed based on scoring them using high-coverage genomes) !!!
		FILE*fixed_suns_file;
		int num_fixed_suns=0;
		char sun_containing_contig_name[51];
		sun_containing_contig_name[50]='\0';
		char note[51];
		note[50]='\0';
		long contig_sun_loc,master_sun_loc;
		double sun_score;
		fixed_suns_file=fopen(*(argv+3+2*num_gene_families+num_prev_gene_family_files+num_contigs_per_family[i]),"r");
		fgetpos(fixed_suns_file,&pos);
		while(fscanf(fixed_suns_file,"%s %ld %ld %lf %s",sun_containing_contig_name,&contig_sun_loc,&master_sun_loc,&sun_score,note)==5)
		{
			num_fixed_suns++;
		}
		fsetpos(fixed_suns_file,&pos);
		long sun_locations[num_fixed_suns];
		k=0;
		while(fscanf(fixed_suns_file,"%s %ld %ld %lf %s",sun_containing_contig_name,&contig_sun_loc,&master_sun_loc,&sun_score,note)==5)
    {
      sun_locations[k]=master_sun_loc;
			k++;
    }
		fclose(fixed_suns_file);

		//input exons
		FILE*exons_file;
		int num_exons=0;
		exons_file=fopen(*(argv+3+2*num_gene_families+num_prev_gene_family_files+num_contigs_per_family[i]+1),"r");
		fgetpos(exons_file,&pos);
		while(fscanf(exons_file,"%s %ld %ld",sun_containing_contig_name,&contig_sun_loc,&master_sun_loc)==3)
		{
			num_exons++;
		}
		fsetpos(exons_file,&pos);
		long exon_starts[num_exons],exon_ends[num_exons];
		k=0;
		while(fscanf(exons_file,"%s %ld %ld",sun_containing_contig_name,&contig_sun_loc,&master_sun_loc)==3)
    {
      exon_starts[k]=contig_sun_loc;
			exon_ends[k]=master_sun_loc;
			k++;
    }
		fclose(exons_file);

		//process mip target information for mips targeting that specific gene family
		long master_start,master_end,contig_start;
		int targets_sun,targets_exon;
		int first_contig_seq;
		int has_variation_between_paralogs;
		char*contig_start_ptr,*contig_start_loc_ptr,*other_arm_loc_ptr,*search_ptr;
		char oldbase,newbase;
		int num_matches;
		char temp[51];
		char contig_target_seqs[num_contigs_per_family[i]][153];
		for(j=0;j<(num_contigs_per_family[i]);j++)
		{
			contig_target_seqs[j][152]='\0';
		}
		for(j=0;j<num_mips_per_family[i];j++)
		{
			master_start=0;
			master_end=0;
			targets_sun=0;
			targets_exon=0;
			fscanf(mip_arm_locations_file,"%s %ld %ld %ld %ld",master_sequence_target,&ext_start,&ext_end,&lig_start,&lig_end);
			
			//store information on target sequence
			strncpy(mip_targets[m].target_sequence,master_sequence_target,50);
			
			//calculate master sequence start and end coordinates of target sequence, mip orientation, and length of first hybridization arm
			if(ext_start<lig_start)
			{
				mip_targets[m].mip_orientation='+';
				master_start=ext_start;
				master_end=lig_end;
				mip_targets[m].first_arm_length=ext_end-ext_start+1;
			}
			else
			{
				mip_targets[m].mip_orientation='-';
				master_start=lig_start;
				master_end=ext_end;
				mip_targets[m].first_arm_length=lig_end-lig_start+1;
			}
			sprintf(mip_targets[m].master_loc,"%ld,%ld",master_start,master_end);
			
			//calculate local coordinates within 152 bp target of any fixed SUNs and determine whether or not target contains any fixed suns
			for(k=0;k<num_fixed_suns;k++)
			{
				if((sun_locations[k]>master_start)&&(sun_locations[k]<master_end))
				{
					targets_sun=1;
					sprintf(temp,"%ld,",sun_locations[k]-master_start+1);
					strcat(mip_targets[m].fixed_sun_locs,temp);
				}
			}
			//if(!(targets_sun)) !!! changed from original version because SUN associated with new MIP designs were not actually scored (scores are "fake" scores) !!!
			//{
			sprintf(mip_targets[m].fixed_sun_locs,"0,\0");
			//}

			//determine whether or not target contains exonic and/or splice site sequence
			for(k=0;k<num_exons;k++)
			{
				if((((exon_starts[k]-2)>master_start)&&((exon_starts[k]-2)<master_end))||(((exon_ends[k]+2)>master_start)&&((exon_ends[k]+2)<master_end))||((master_start>exon_starts[k]-2)&&(master_start<exon_ends[k]+2))) // +/- 2 to account for splice site sequences
				{
					targets_exon=1;
				}
			}

			//use information about SUN and exonic sequence content of target to specify mip type (SUN-targeting, exon-targeting, or both)
			if((targets_sun)&&(targets_exon))
			{
				mip_targets[m].mip_type='B';
			}
			else if(targets_sun)
			{
				mip_targets[m].mip_type='S';
			}
			else if(targets_exon)
			{
				mip_targets[m].mip_type='E';
			}
			else
			{
				mip_targets[m].mip_type='N';
				printf("WARNING! MIP targeting %s bases %ld-%ld targets neither a SUN nor an exon!\n",master_sequence_target,master_start,master_end);
			}

			//get first hybridization arm sequence and search contig sequences for it to calculate contig start coordinates of target sequence; read target sequences from contigs into character arrays
			// !!! differs from "detail_mip_targets.c" as it looks for both hybridization arm sequences !!!
			char first_arm_seq[mip_targets[m].first_arm_length+1];
			char other_arm_seq[40-mip_targets[m].first_arm_length+1];
			first_arm_seq[mip_targets[m].first_arm_length]='\0';
			other_arm_seq[40-mip_targets[m].first_arm_length]='\0';
			strncpy(first_arm_seq,master_sequence+master_start-1,mip_targets[m].first_arm_length);
			strncpy(other_arm_seq,master_sequence+master_start-1+mip_targets[m].first_arm_length+112,40-mip_targets[m].first_arm_length);
			int target_found;
			for(k=0;k<num_contigs_per_family[i];k++)
			{
				target_found=0;
				contig_start_ptr=&(contig_sequences[k][0]);
				search_ptr=contig_start_ptr;
				while((contig_start_loc_ptr=strstr(search_ptr,first_arm_seq))!=NULL)
				{
					search_ptr=contig_start_loc_ptr+1;
					other_arm_loc_ptr=strstr(search_ptr,other_arm_seq);
					if((other_arm_loc_ptr!=NULL)&&(other_arm_loc_ptr-(contig_start_loc_ptr+mip_targets[m].first_arm_length)==112))
					{
						target_found=1;
						break;
					}
				}
				if(!(target_found))
				{
					sprintf(temp,"%s:0;",contig_sequence_names[k]);
					strcat(mip_targets[m].contig_start_locs,temp);
					for(n=0;n<152;n++)
					{
						contig_target_seqs[k][n]='-';
					}
				}
				else
				{
					sprintf(temp,"%s:%ld;",contig_sequence_names[k],(long)(contig_start_loc_ptr-contig_start_ptr+1));
					strcat(mip_targets[m].contig_start_locs,temp);
					strncpy(contig_target_seqs[k],contig_start_loc_ptr,152);
				}
			}

			//calculate local coordinates within 152 bp target of any sites where variation in the form of substitution between paralogs is observed
			has_variation_between_paralogs=0;
			for(n=0;n<152;n++)
			{
				first_contig_seq=0;
				oldbase=contig_target_seqs[first_contig_seq][n];
				while(oldbase=='-')
				{
					first_contig_seq++;
					oldbase=contig_target_seqs[first_contig_seq][n];
				}
				for(k=1;k<num_contigs_per_family[i];k++)
				{
					newbase=contig_target_seqs[k][n];
					if((newbase!=oldbase)&&(newbase!='-'))
					{
						sprintf(temp,"%d,",n+1);
						strcat(mip_targets[m].variant_locs,temp);
						has_variation_between_paralogs=1;
						break;
					}
				}
			}
			if(!(has_variation_between_paralogs))
			{
				sprintf(mip_targets[m].variant_locs,"0,");
			}

			//determine whether each contig sequence confers paralog-specificity and output this information in the form of a string of letters
			char seqtype='A';
			char specvec[num_contigs_per_family[i]+1];
			for(k=0;k<=num_contigs_per_family[i];k++)
			{
				specvec[k]='\0';
			}
			for(k=0;k<num_contigs_per_family[i];k++)
			{
				if(specvec[k]>='A')
				{
					continue;
				}
				first_contig_seq=k;
				for(n=0;n<num_contigs_per_family[i];n++)
				{
					if(strcmp(contig_target_seqs[first_contig_seq],contig_target_seqs[n])==0)
					{
						specvec[n]=seqtype;
					}
				}
				seqtype++;
			}
			strcat(mip_targets[m].specificity,specvec);
			m++;
		}
	}
	
	//print information on each mip target to output file
	FILE*out;
	char extension[12]=".miptargets";
	char output_file_name[62];
	sprintf(output_file_name,"%s",*(argv+(argc-1)));
	strcat(output_file_name,extension);
	out=fopen(output_file_name,"w");
	for(m=0;m<num_mips;m++)
	{
		fprintf(out,"%s\t%s\t%s\t%c\t%s\t%s\t%s\t%c\t%d\n",mip_targets[m].target_sequence,mip_targets[m].master_loc,mip_targets[m].contig_start_locs,mip_targets[m].mip_type,mip_targets[m].variant_locs,mip_targets[m].fixed_sun_locs,mip_targets[m].specificity,mip_targets[m].mip_orientation,mip_targets[m].first_arm_length);
	}

	//print any targeting conflicts to standard output
	long start,start2,end,end2;
	for(m=0;m<num_mips;m++)
	{
		start=strtol(mip_targets[m].master_loc,NULL,10);
		end=strtol(strchr(mip_targets[m].master_loc,',')+1,NULL,10);
		for(n=0;n<num_mips;n++)
		{
			start2=strtol(mip_targets[n].master_loc,NULL,10);
    	end2=strtol(strchr(mip_targets[n].master_loc,',')+1,NULL,10);
			if((strcmp(mip_targets[m].target_sequence,mip_targets[n].target_sequence)==0)&&(m!=n)&&(mip_targets[m].mip_orientation==mip_targets[n].mip_orientation)&&(((start>=start2)&&(start<=end2))||((end>=start2)&&(end<=end2))))
			{
				printf("Targeting conflict! MIP targeting %s:%ld-%ld conflicts with MIP targeting %s:%ld-%ld\n",mip_targets[m].target_sequence,start,end,mip_targets[n].target_sequence,start2,end2);
			}
		}
	}

	//clean up and exit
	fclose(mip_arm_locations_file);
	fclose(out);
	return 0;
}
